#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <vector>

// Implements the small-SNN kernel-width strategy used by Seurat 5.5.0 WNN
// without depending on Seurat's private compiled interface. Ordering by shared
// neighbor count is equivalent to ordering by its Jaccard transformation.
extern "C" SEXP multiomeR_wnn_small_snn_bandwidth(
    SEXP embedding_sexp,
    SEXP knn_idx_sexp,
    SEXP k_sexp,
    SEXP nearest_dist_sexp) {
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
          embedding[cell + cell_count * dimension] -
          embedding[entry.first + cell_count * dimension];
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

    if ((cell & 0x3ff) == 0) {
      R_CheckUserInterrupt();
    }
  }

  UNPROTECT(1);
  return bandwidth_sexp;
}
