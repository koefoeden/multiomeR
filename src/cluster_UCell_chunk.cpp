#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <functional>
#include <vector>

// Cluster sums of capped UCell rank gaps for one sparse count chunk.
// Zero counts share the average rank (rows + nonzero + 1) / 2, so each cell is
// ranked from its positive entries only. Means use a long double sum divided by
// the gene count, as R's colMeans(); cluster sums add cells in chunk order in
// double, as R's rowsum(). Adding exact zeros is skipped because it is an
// identity, so the results equal the dense R implementation bit for bit.
namespace {

struct Component {
  const int *index = nullptr;  // genes x (controls + 1), column-major, 0-based
  int genes = 0;
  double correction = 1.0;
};

inline double component_mean(const Component &component, const double *gap, int replicate) {
  if (!component.genes) return 0.0;
  const int *index = component.index + static_cast<R_xlen_t>(replicate) * component.genes;
  long double sum = 0.0;
  for (int gene = 0; gene < component.genes; ++gene) sum += gap[index[gene]];
  if (sum == 0.0) return 0.0;
  sum /= component.genes;
  return static_cast<double>(sum) / component.correction;
}

Component read_component(SEXP index, SEXP correction, int replicates, int genes) {
  Component component;
  if (Rf_isNull(index)) return component;
  if (TYPEOF(index) != INTSXP || Rf_length(index) % replicates) Rf_error("Invalid UCell component index");
  component.index = INTEGER(index);
  component.genes = Rf_length(index) / replicates;
  component.correction = Rf_asReal(correction);
  for (R_xlen_t k = 0; k < Rf_xlength(index); ++k) {
    if (component.index[k] < 0 || component.index[k] >= genes) Rf_error("UCell component index out of range");
  }
  return component;
}

}  // namespace

extern "C" SEXP multiomeR_UCell_chunk(SEXP p_sexp, SEXP i_sexp, SEXP x_sexp, SEXP rows_sexp,
                                      SEXP gene_position_sexp, SEXP genes_sexp, SEXP cluster_sexp,
                                      SEXP clusters_sexp, SEXP labels_sexp, SEXP replicates_sexp,
                                      SEXP max_rank_sexp) {
  const int *p = INTEGER(p_sexp), *row = INTEGER(i_sexp), *gene_position = INTEGER(gene_position_sexp);
  const int *cluster = INTEGER(cluster_sexp);
  const double *x = REAL(x_sexp);
  const int rows = Rf_asInteger(rows_sexp), genes = Rf_asInteger(genes_sexp);
  const int clusters = Rf_asInteger(clusters_sexp), replicates = Rf_asInteger(replicates_sexp);
  const int cells = Rf_length(cluster_sexp), labels = Rf_length(labels_sexp);
  const double max_rank = Rf_asReal(max_rank_sexp);
  if (Rf_length(p_sexp) != cells + 1 || Rf_length(gene_position_sexp) != rows) Rf_error("Invalid UCell chunk dimensions");

  std::vector<Component> positive(labels), negative(labels);
  std::vector<int> is_signed(labels);
  for (int label = 0; label < labels; ++label) {
    SEXP spec = VECTOR_ELT(labels_sexp, label);
    positive[label] = read_component(VECTOR_ELT(spec, 0), VECTOR_ELT(spec, 1), replicates, genes);
    negative[label] = read_component(VECTOR_ELT(spec, 2), VECTOR_ELT(spec, 3), replicates, genes);
    is_signed[label] = negative[label].genes > 0;
  }

  SEXP rank_sum = PROTECT(Rf_allocMatrix(REALSXP, genes, clusters));
  SEXP detected_sum = PROTECT(Rf_allocMatrix(REALSXP, genes, clusters));
  SEXP signed_sum = PROTECT(Rf_allocVector(VECSXP, labels));
  SEXP cell_scores = PROTECT(Rf_allocMatrix(REALSXP, cells, labels));
  double *rank_out = REAL(rank_sum), *detected_out = REAL(detected_sum), *cell_out = REAL(cell_scores);
  std::fill(rank_out, rank_out + static_cast<R_xlen_t>(genes) * clusters, 0.0);
  std::fill(detected_out, detected_out + static_cast<R_xlen_t>(genes) * clusters, 0.0);
  std::vector<double *> signed_out(labels, nullptr);
  for (int label = 0; label < labels; ++label) {
    if (!is_signed[label]) continue;
    SET_VECTOR_ELT(signed_sum, label, Rf_allocMatrix(REALSXP, replicates, clusters));
    signed_out[label] = REAL(VECTOR_ELT(signed_sum, label));
    std::fill(signed_out[label], signed_out[label] + static_cast<R_xlen_t>(replicates) * clusters, 0.0);
  }

  std::vector<double> gap(genes), values;
  for (int cell = 0; cell < cells; ++cell) {
    const int start = p[cell], end = p[cell + 1];
    const int group = cluster[cell];
    if (group < 0 || group >= clusters) Rf_error("UCell cluster index out of range");
    values.clear();
    for (int entry = start; entry < end; ++entry) {
      if (!(x[entry] >= 0.0)) Rf_error("UCell counts must be non-negative and finite");
      if (x[entry] > 0.0) values.push_back(x[entry]);
    }
    std::sort(values.begin(), values.end(), std::greater<double>());
    const double zero_rank = (rows + static_cast<double>(values.size()) + 1.0) / 2.0;
    std::fill(gap.begin(), gap.end(), (max_rank - std::min(zero_rank, max_rank)) / max_rank);
    double *detected_group = detected_out + static_cast<R_xlen_t>(group) * genes;
    for (int entry = start; entry < end; ++entry) {
      const int position = gene_position[row[entry]];
      if (position < 0 || x[entry] == 0.0) continue;
      const auto first = std::lower_bound(values.begin(), values.end(), x[entry], std::greater<double>());
      const auto last = std::upper_bound(first, values.end(), x[entry], std::greater<double>());
      const double rank = ((first - values.begin()) + 1.0 + (last - values.begin())) / 2.0;
      gap[position] = (max_rank - std::min(rank, max_rank)) / max_rank;
      detected_group[position] += 1.0;
    }
    double *rank_group = rank_out + static_cast<R_xlen_t>(group) * genes;
    for (int gene = 0; gene < genes; ++gene) {
      if (gap[gene] != 0.0) rank_group[gene] += gap[gene];
    }
    for (int label = 0; label < labels; ++label) {
      if (!is_signed[label]) {
        cell_out[cell + static_cast<R_xlen_t>(label) * cells] = component_mean(positive[label], gap.data(), 0);
        continue;
      }
      double *signed_group = signed_out[label] + static_cast<R_xlen_t>(group) * replicates;
      for (int replicate = 0; replicate < replicates; ++replicate) {
        double score = component_mean(positive[label], gap.data(), replicate) -
                       component_mean(negative[label], gap.data(), replicate);
        if (!(score > 0.0)) score = 0.0;  // pmax(0, score)
        if (score != 0.0) signed_group[replicate] += score;
        if (!replicate) cell_out[cell + static_cast<R_xlen_t>(label) * cells] = score;
      }
    }
  }

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(result, 0, rank_sum);
  SET_VECTOR_ELT(result, 1, detected_sum);
  SET_VECTOR_ELT(result, 2, signed_sum);
  SET_VECTOR_ELT(result, 3, cell_scores);
  UNPROTECT(5);
  return result;
}
