#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <vector>

extern "C" SEXP multiomeR_scavenge_permutation_statistics(
    SEXP column_ptr_sexp,
    SEXP row_idx_sexp,
    SEXP values_sexp,
    SEXP seed_offsets_sexp,
    SEXP seed_indices_sexp,
    SEXP chunk_offsets_sexp,
    SEXP observed_scores_sexp,
    SEXP cluster_offsets_sexp,
    SEXP cluster_cell_indices_sexp,
    SEXP restart_prob_sexp,
    SEXP stationary_cutoff_sexp,
    SEXP max_iter_sexp,
    SEXP cores_sexp) {
  if (TYPEOF(column_ptr_sexp) != INTSXP ||
      TYPEOF(row_idx_sexp) != INTSXP ||
      TYPEOF(values_sexp) != REALSXP ||
      TYPEOF(seed_offsets_sexp) != INTSXP ||
      TYPEOF(seed_indices_sexp) != INTSXP ||
      TYPEOF(chunk_offsets_sexp) != INTSXP ||
      TYPEOF(observed_scores_sexp) != REALSXP ||
      TYPEOF(cluster_offsets_sexp) != INTSXP ||
      TYPEOF(cluster_cell_indices_sexp) != INTSXP ||
      TYPEOF(restart_prob_sexp) != REALSXP ||
      TYPEOF(stationary_cutoff_sexp) != REALSXP ||
      TYPEOF(max_iter_sexp) != INTSXP ||
      TYPEOF(cores_sexp) != INTSXP ||
      Rf_length(restart_prob_sexp) != 1 ||
      Rf_length(stationary_cutoff_sexp) != 1 ||
      Rf_length(max_iter_sexp) != 1 ||
      Rf_length(cores_sexp) != 1) {
    Rf_error("Invalid SCAVENGE native-kernel input type");
  }

  const int cell_count = Rf_length(observed_scores_sexp);
  const int edge_count = Rf_length(row_idx_sexp);
  const int permutation_count = Rf_length(seed_offsets_sexp) - 1;
  const int chunk_count = Rf_length(chunk_offsets_sexp) - 1;
  const int cluster_count = Rf_length(cluster_offsets_sexp) - 1;
  const double restart_prob = REAL(restart_prob_sexp)[0];
  const double stationary_cutoff = REAL(stationary_cutoff_sexp)[0];
  const int max_iter = INTEGER(max_iter_sexp)[0];
  const int cores = INTEGER(cores_sexp)[0];

  if (Rf_length(column_ptr_sexp) != cell_count + 1 ||
      Rf_length(values_sexp) != edge_count ||
      permutation_count < 1 ||
      chunk_count < 1 ||
      cluster_count < 1 ||
      INTEGER(column_ptr_sexp)[0] != 0 ||
      INTEGER(column_ptr_sexp)[cell_count] != edge_count ||
      INTEGER(seed_offsets_sexp)[0] != 0 ||
      INTEGER(seed_offsets_sexp)[permutation_count] !=
          Rf_length(seed_indices_sexp) ||
      INTEGER(chunk_offsets_sexp)[0] != 0 ||
      INTEGER(chunk_offsets_sexp)[chunk_count] != permutation_count ||
      INTEGER(cluster_offsets_sexp)[0] != 0 ||
      INTEGER(cluster_offsets_sexp)[cluster_count] !=
          Rf_length(cluster_cell_indices_sexp) ||
      !(restart_prob > 0.0 && restart_prob < 1.0) ||
      stationary_cutoff < 0.0 ||
      max_iter < 1 ||
      cores < 1) {
    Rf_error("Invalid SCAVENGE native-kernel input value");
  }
  for (int permutation = 0; permutation < permutation_count; ++permutation) {
    if (INTEGER(seed_offsets_sexp)[permutation] >=
        INTEGER(seed_offsets_sexp)[permutation + 1]) {
      Rf_error("Every SCAVENGE permutation must contain at least one seed");
    }
  }
  for (int chunk = 0; chunk < chunk_count; ++chunk) {
    if (INTEGER(chunk_offsets_sexp)[chunk] >=
        INTEGER(chunk_offsets_sexp)[chunk + 1]) {
      Rf_error("Every SCAVENGE work chunk must contain a permutation");
    }
  }
  int max_cluster_size = 0;
  for (int cluster = 0; cluster < cluster_count; ++cluster) {
    const int cluster_size =
        INTEGER(cluster_offsets_sexp)[cluster + 1] -
        INTEGER(cluster_offsets_sexp)[cluster];
    if (cluster_size < 1) {
      Rf_error("Every SCAVENGE cluster must contain at least one cell");
    }
    max_cluster_size = std::max(max_cluster_size, cluster_size);
  }
  for (int i = 0; i < Rf_length(seed_indices_sexp); ++i) {
    const int cell = INTEGER(seed_indices_sexp)[i];
    if (cell < 0 || cell >= cell_count) {
      Rf_error("SCAVENGE seed index is out of bounds");
    }
  }
  for (int i = 0; i < Rf_length(cluster_cell_indices_sexp); ++i) {
    const int cell = INTEGER(cluster_cell_indices_sexp)[i];
    if (cell < 0 || cell >= cell_count) {
      Rf_error("SCAVENGE cluster cell index is out of bounds");
    }
  }

  const int *column_ptr = INTEGER(column_ptr_sexp);
  const int *row_idx = INTEGER(row_idx_sexp);
  const double *values = REAL(values_sexp);
  const int *seed_offsets = INTEGER(seed_offsets_sexp);
  const int *seed_indices = INTEGER(seed_indices_sexp);
  const int *chunk_offsets = INTEGER(chunk_offsets_sexp);
  const double *observed_scores = REAL(observed_scores_sexp);
  const int *cluster_offsets = INTEGER(cluster_offsets_sexp);
  const int *cluster_cell_indices =
      INTEGER(cluster_cell_indices_sexp);

  SEXP exceedance_counts = PROTECT(
      Rf_allocMatrix(INTSXP, cell_count, chunk_count));
  std::fill(
      INTEGER(exceedance_counts),
      INTEGER(exceedance_counts) +
          static_cast<R_xlen_t>(cell_count) * chunk_count,
      0);
  int *counts = INTEGER(exceedance_counts);
  SEXP cluster_statistics = PROTECT(
      Rf_allocMatrix(REALSXP, cluster_count, permutation_count));
  double *cluster_output = REAL(cluster_statistics);
  const double walk_prob = 1.0 - restart_prob;

#ifdef _OPENMP
#pragma omp parallel num_threads(cores)
#endif
  {
    std::vector<double> score(cell_count);
    std::vector<double> next_score(cell_count);
    std::vector<double> cluster_scores(max_cluster_size);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int chunk = 0; chunk < chunk_count; ++chunk) {
      int *chunk_counts =
          counts + static_cast<R_xlen_t>(chunk) * cell_count;

      for (
          int permutation = chunk_offsets[chunk];
          permutation < chunk_offsets[chunk + 1];
          ++permutation) {
        std::fill(score.begin(), score.end(), 0.0);
        const int seed_begin = seed_offsets[permutation];
        const int seed_end = seed_offsets[permutation + 1];
        const double seed_weight =
            1.0 / static_cast<double>(seed_end - seed_begin);
        for (int seed = seed_begin; seed < seed_end; ++seed) {
          score[seed_indices[seed]] = seed_weight;
        }

        for (int iteration = 0; iteration < max_iter; ++iteration) {
          std::fill(next_score.begin(), next_score.end(), 0.0);
          for (int column = 0; column < cell_count; ++column) {
            const double column_score = score[column];
            if (column_score == 0.0) {
              continue;
            }
            for (
                int edge = column_ptr[column];
                edge < column_ptr[column + 1];
                ++edge) {
              next_score[row_idx[edge]] += values[edge] * column_score;
            }
          }

          long double delta = 0.0L;
          for (int cell = 0; cell < cell_count; ++cell) {
            next_score[cell] *= walk_prob;
          }
          for (int seed = seed_begin; seed < seed_end; ++seed) {
            next_score[seed_indices[seed]] += restart_prob * seed_weight;
          }
          for (int cell = 0; cell < cell_count; ++cell) {
            delta += std::abs(next_score[cell] - score[cell]);
          }
          score.swap(next_score);
          if (delta <= stationary_cutoff) {
            break;
          }
        }

        for (int cell = 0; cell < cell_count; ++cell) {
          chunk_counts[cell] +=
              static_cast<int>(score[cell] > observed_scores[cell]);
        }
        for (int cluster = 0; cluster < cluster_count; ++cluster) {
          const int cluster_begin = cluster_offsets[cluster];
          const int cluster_size =
              cluster_offsets[cluster + 1] - cluster_begin;
          for (int i = 0; i < cluster_size; ++i) {
            cluster_scores[i] =
                score[cluster_cell_indices[cluster_begin + i]];
          }

          double *scores_begin = cluster_scores.data();
          double *scores_mid = scores_begin + cluster_size / 2;
          std::nth_element(
              scores_begin,
              scores_mid,
              scores_begin + cluster_size);
          double median = *scores_mid;
          if (cluster_size % 2 == 0) {
            median =
                (*std::max_element(scores_begin, scores_mid) + median) /
                2.0;
          }
          cluster_output[
              cluster +
              static_cast<R_xlen_t>(cluster_count) * permutation] =
              median;
        }
      }
    }
  }

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, exceedance_counts);
  SET_VECTOR_ELT(result, 1, cluster_statistics);
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(
      names,
      0,
      Rf_mkChar("cell_exceedance_counts"));
  SET_STRING_ELT(
      names,
      1,
      Rf_mkChar("cluster_statistics"));
  Rf_setAttrib(result, R_NamesSymbol, names);

  UNPROTECT(4);
  return result;
}
