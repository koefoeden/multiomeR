#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Implements the small-SNN kernel-width strategy used by Seurat 5.5.0 WNN
// without depending on Seurat's private compiled interface. Ordering by shared
// neighbor count is equivalent to ordering by its Jaccard transformation.
extern "C" SEXP multiomeR_wnn_small_snn_bandwidth(
    SEXP embedding_sexp,
    SEXP knn_idx_sexp,
    SEXP k_sexp,
    SEXP nearest_dist_sexp,
    SEXP threads_sexp) {
  const int threads = Rf_asInteger(threads_sexp);
  if (TYPEOF(embedding_sexp) != REALSXP ||
      TYPEOF(knn_idx_sexp) != INTSXP ||
      TYPEOF(k_sexp) != INTSXP ||
      TYPEOF(nearest_dist_sexp) != REALSXP ||
      Rf_length(k_sexp) != 1) {
    Rf_error("Invalid WNN bandwidth input type");
  }

  SEXP embedding_dims = Rf_getAttrib(embedding_sexp, R_DimSymbol);
  SEXP knn_dims = Rf_getAttrib(knn_idx_sexp, R_DimSymbol);
  if (Rf_length(embedding_dims) != 2 || Rf_length(knn_dims) != 2) {
    Rf_error("WNN bandwidth inputs must be matrices");
  }

  const int cell_count = INTEGER(embedding_dims)[0];
  const int dimension_count = INTEGER(embedding_dims)[1];
  const int knn_cell_count = INTEGER(knn_dims)[0];
  const int neighbor_count = INTEGER(knn_dims)[1];
  const int k = INTEGER(k_sexp)[0];
  if (cell_count < 2 || dimension_count < 1 ||
      knn_cell_count != cell_count ||
      Rf_length(nearest_dist_sexp) != cell_count ||
      k < 1 || k > neighbor_count) {
    Rf_error("Invalid WNN bandwidth input dimensions");
  }

  const double *embedding = REAL(embedding_sexp);
  // Row-major copy so that each cell's coordinates are contiguous.
  std::vector<double> embedding_t(static_cast<std::size_t>(cell_count) * dimension_count);
  for (int dimension = 0; dimension < dimension_count; ++dimension) {
    for (int cell = 0; cell < cell_count; ++cell) {
      embedding_t[dimension + dimension_count * static_cast<std::size_t>(cell)] =
        embedding[cell + static_cast<std::size_t>(cell_count) * dimension];
    }
  }
  const int *knn_idx = INTEGER(knn_idx_sexp);
  const double *nearest_dist = REAL(nearest_dist_sexp);

  std::vector<std::vector<int>> cells_by_neighbor(cell_count);
  for (int cell = 0; cell < cell_count; ++cell) {
    for (int rank = 0; rank < k; ++rank) {
      const int neighbor = knn_idx[cell + cell_count * rank] - 1;
      if (neighbor < 0 || neighbor >= cell_count) {
        Rf_error("WNN neighbor index is out of bounds");
      }
      cells_by_neighbor[neighbor].push_back(cell);
    }
  }

  SEXP bandwidth_sexp = PROTECT(Rf_allocVector(REALSXP, cell_count));
  double *bandwidth = REAL(bandwidth_sexp);

#pragma omp parallel for num_threads(threads) schedule(dynamic, 256)
  for (int cell = 0; cell < cell_count; ++cell) {
    std::unordered_map<int, int> shared_neighbor_counts;
    shared_neighbor_counts.reserve(static_cast<std::size_t>(k) * k);
    for (int rank = 0; rank < k; ++rank) {
      const int neighbor = knn_idx[cell + cell_count * rank] - 1;
      for (const int other_cell : cells_by_neighbor[neighbor]) {
        ++shared_neighbor_counts[other_cell];
      }
    }

    const int selected_count = std::min(
      k,
      static_cast<int>(shared_neighbor_counts.size())
    );
    std::vector<int> counts;
    counts.reserve(shared_neighbor_counts.size());
    for (const auto &entry : shared_neighbor_counts) {
      counts.push_back(entry.second);
    }
    std::nth_element(
      counts.begin(),
      counts.begin() + selected_count - 1,
      counts.end()
    );
    const int count_threshold = counts[selected_count - 1];

    std::vector<double> distances;
    distances.reserve(shared_neighbor_counts.size());
    for (const auto &entry : shared_neighbor_counts) {
      if (entry.second > count_threshold) {
        continue;
      }

      double squared_distance = 0.0;
      for (int dimension = 0; dimension < dimension_count; ++dimension) {
        const double delta =
          embedding_t[dimension + dimension_count * static_cast<std::size_t>(cell)] -
          embedding_t[dimension + dimension_count * static_cast<std::size_t>(entry.first)];
        squared_distance += delta * delta;
      }
      distances.push_back(std::max(
        std::sqrt(squared_distance) - nearest_dist[cell],
        0.0
      ));
    }

    if (static_cast<int>(distances.size()) > selected_count) {
      std::sort(distances.rbegin(), distances.rend());
      distances.resize(selected_count);
    }
    bandwidth[cell] = std::accumulate(
      distances.begin(),
      distances.end(),
      0.0
    ) / distances.size();

  }

  UNPROTECT(1);
  return bandwidth_sexp;
}

// Weighted WNN neighbour selection, reproducing the per-cell R reference
// bit for bit: candidates are the first occurrences of each modality's
// neighbours 2..candidate_k+1 in modality order; each Euclidean distance
// accumulates squared double differences in long double over dimensions in
// order, as rowSums() does; scores add modalities left to right; the top k
// follow a stable decreasing order, as order(decreasing = TRUE) does.
// Embeddings arrive transposed (dimensions x cells) so that each cell's
// coordinates are contiguous; values and summation order are unchanged.
extern "C" SEXP multiomeR_wnn_select(
    SEXP embeddings_sexp,
    SEXP knn_idx_sexp,
    SEXP nearest_dist_sexp,
    SEXP sigma_sexp,
    SEXP weight_sexp,
    SEXP candidate_k_sexp,
    SEXP k_sexp,
    SEXP threads_sexp) {
  const int modality_count = Rf_length(embeddings_sexp);
  const int candidate_k = Rf_asInteger(candidate_k_sexp);
  const int k = Rf_asInteger(k_sexp);
  const int threads = Rf_asInteger(threads_sexp);
  const int cell_count = Rf_ncols(VECTOR_ELT(embeddings_sexp, 0));

  std::vector<const double *> embedding(modality_count), nearest_dist(modality_count),
      sigma(modality_count), weight(modality_count);
  std::vector<const int *> knn_idx(modality_count);
  std::vector<int> dimension_count(modality_count);
  for (int m = 0; m < modality_count; ++m) {
    SEXP e = VECTOR_ELT(embeddings_sexp, m), idx = VECTOR_ELT(knn_idx_sexp, m);
    if (TYPEOF(e) != REALSXP || TYPEOF(idx) != INTSXP || Rf_ncols(e) != cell_count ||
        Rf_nrows(idx) != cell_count || Rf_ncols(idx) < candidate_k + 1 ||
        Rf_length(VECTOR_ELT(nearest_dist_sexp, m)) != cell_count ||
        Rf_length(VECTOR_ELT(sigma_sexp, m)) != cell_count ||
        Rf_length(VECTOR_ELT(weight_sexp, m)) != cell_count) {
      Rf_error("Invalid WNN selection input");
    }
    embedding[m] = REAL(e);
    dimension_count[m] = Rf_nrows(e);
    knn_idx[m] = INTEGER(idx);
    nearest_dist[m] = REAL(VECTOR_ELT(nearest_dist_sexp, m));
    sigma[m] = REAL(VECTOR_ELT(sigma_sexp, m));
    weight[m] = REAL(VECTOR_ELT(weight_sexp, m));
  }
  if (k < 1 || k > candidate_k) {
    Rf_error("Invalid WNN selection k");
  }

  SEXP nn_idx_sexp = PROTECT(Rf_allocMatrix(INTSXP, cell_count, k));
  SEXP nn_dist_sexp = PROTECT(Rf_allocMatrix(REALSXP, cell_count, k));
  int *nn_idx = INTEGER(nn_idx_sexp);
  double *nn_dist = REAL(nn_dist_sexp);
  const R_xlen_t n = cell_count;
  int invalid = 0;

#pragma omp parallel num_threads(threads)
  {
    std::vector<int> seen(cell_count, -1);
    std::vector<int> candidates;
    std::vector<double> score;
    std::vector<int> order;
    candidates.reserve(static_cast<std::size_t>(modality_count) * candidate_k);

#pragma omp for schedule(dynamic, 256)
    for (int cell = 0; cell < cell_count; ++cell) {
      candidates.clear();
      for (int m = 0; m < modality_count; ++m) {
        for (int rank = 1; rank <= candidate_k; ++rank) {
          const int neighbor = knn_idx[m][cell + n * rank] - 1;
          if (neighbor < 0 || neighbor >= cell_count) {
#pragma omp atomic write
            invalid = 1;
            continue;
          }
          if (seen[neighbor] != cell) {
            seen[neighbor] = cell;
            candidates.push_back(neighbor);
          }
        }
      }

      const int candidate_count = static_cast<int>(candidates.size());
      score.assign(candidate_count, 0.0);
      for (int m = 0; m < modality_count; ++m) {
        const int dims = dimension_count[m];
        const double *query = embedding[m] + static_cast<R_xlen_t>(dims) * cell;
        for (int c = 0; c < candidate_count; ++c) {
          const double *x = embedding[m] + static_cast<R_xlen_t>(dims) * candidates[c];
          long double squared_distance = 0.0L;
          for (int dim = 0; dim < dims; ++dim) {
            const double delta = query[dim] - x[dim];
            const double squared = delta * delta;
            squared_distance += squared;
          }
          double distance = std::sqrt(static_cast<double>(squared_distance)) - nearest_dist[m][cell];
          if (0.0 > distance) distance = 0.0;  // pmax(distance, 0)
          const double term = std::exp(-1.0 * (distance / sigma[m][cell])) * weight[m][cell];
          score[c] = m == 0 ? term : score[c] + term;
        }
      }

      order.resize(candidate_count);
      std::iota(order.begin(), order.end(), 0);
      std::partial_sort(order.begin(), order.begin() + k, order.end(), [&](int a, int b) {
        return score[a] > score[b] || (score[a] == score[b] && a < b);
      });
      for (int rank = 0; rank < k; ++rank) {
        const double s = score[order[rank]];
        double half_gap = (1 - s) / 2;
        if (0.0 > half_gap) half_gap = 0.0;  // pmax(., 0)
        if (1.0 < half_gap) half_gap = 1.0;  // pmin(., 1)
        nn_idx[cell + n * rank] = candidates[order[rank]] + 1;
        nn_dist[cell + n * rank] = std::sqrt(half_gap);
      }
    }
  }
  if (invalid) {
    UNPROTECT(2);
    Rf_error("WNN neighbor index is out of bounds");
  }

  SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, nn_idx_sexp);
  SET_VECTOR_ELT(result, 1, nn_dist_sexp);
  UNPROTECT(3);
  return result;
}
