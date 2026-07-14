#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <climits>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

// BPCells does not export its FragmentLoader header. This interface mirrors the
// pinned BPCells revision checked by the R wrapper so the pipeline can stream
// fragments without materializing them as GRanges. The mirrored interface is
// from BPCells (Copyright 2021 BPCells contributors), used under its MIT terms.
class FragmentLoaderABI {
 public:
  virtual ~FragmentLoaderABI() = default;
  virtual bool isSeekable() const = 0;
  virtual void seek(uint32_t chr_id, uint32_t base) = 0;
  virtual void restart() = 0;
  virtual int chrCount() const = 0;
  virtual int cellCount() const = 0;
  virtual const char *chrNames(uint32_t chr_id) = 0;
  virtual const char *cellNames(uint32_t cell_id) = 0;
  virtual bool nextChr() = 0;
  virtual uint32_t currentChr() const = 0;
  virtual bool load() = 0;
  virtual uint32_t capacity() const = 0;
  virtual uint32_t *cellData() = 0;
  virtual uint32_t *startData() = 0;
  virtual uint32_t *endData() = 0;
};

struct Locus {
  uint32_t chr;
  uint32_t start;
  uint32_t end;
  uint32_t cell;
};

static FragmentLoaderABI *get_loader(SEXP external_pointer) {
  if (TYPEOF(external_pointer) != EXTPTRSXP) {
    Rf_error("Expected a BPCells fragments external pointer");
  }
  auto holder = static_cast<std::unique_ptr<FragmentLoaderABI> *>(
      R_ExternalPtrAddr(external_pointer));
  if (holder == nullptr || holder->get() == nullptr) {
    Rf_error("BPCells fragments external pointer is NULL");
  }
  return holder->get();
}

extern "C" SEXP multiomeR_bpcells_fragment_counts(
    SEXP external_pointer, SEXP n_cells_sexp) {
  if (TYPEOF(n_cells_sexp) != INTSXP || Rf_xlength(n_cells_sexp) != 1) {
    Rf_error("`n_cells` must be one integer");
  }

  FragmentLoaderABI *loader = get_loader(external_pointer);
  int n_cells = INTEGER(n_cells_sexp)[0];
  if (n_cells < 0) {
    Rf_error("`n_cells` must be non-negative");
  }

  SEXP output = PROTECT(Rf_allocVector(INTSXP, n_cells));
  int *counts = INTEGER(output);
  std::fill(counts, counts + n_cells, 0);

  try {
    loader->restart();
    uint64_t fragments_seen = 0;
    while (loader->nextChr()) {
      while (loader->load()) {
        uint32_t capacity = loader->capacity();
        uint32_t *cell = loader->cellData();
        for (uint32_t idx = 0; idx < capacity; ++idx) {
          if (cell[idx] >= static_cast<uint32_t>(n_cells)) {
            throw std::runtime_error("BPCells fragment cell index is out of range");
          }
          if (counts[cell[idx]] < INT_MAX) {
            counts[cell[idx]] += 1;
          }
        }
        fragments_seen += capacity;
        if ((fragments_seen & 0x3ffff) == 0) {
          R_CheckUserInterrupt();
        }
      }
    }
  } catch (const std::exception &exception) {
    UNPROTECT(1);
    Rf_error("%s", exception.what());
  }

  UNPROTECT(1);
  return output;
}

static void append_coverage_loci(
    const std::vector<std::pair<uint32_t, uint32_t>> &fragments,
    uint32_t chromosome,
    uint32_t cell,
    std::vector<Locus> &loci) {
  if (fragments.empty()) {
    return;
  }

  std::vector<std::pair<uint32_t, int>> events;
  events.reserve(fragments.size() * 2);
  for (const auto &fragment : fragments) {
    if (fragment.second <= fragment.first) {
      continue;
    }
    events.push_back(std::make_pair(fragment.first, 1));
    events.push_back(std::make_pair(fragment.second, -1));
  }
  if (events.empty()) {
    return;
  }

  std::sort(
      events.begin(),
      events.end(),
      [](const std::pair<uint32_t, int> &left,
         const std::pair<uint32_t, int> &right) {
        if (left.first != right.first) {
          return left.first < right.first;
        }
        return left.second < right.second;
      });

  uint32_t previous_position = events.front().first;
  uint32_t run_start = 0;
  int coverage = 0;
  bool in_run = false;
  size_t event_idx = 0;

  while (event_idx < events.size()) {
    uint32_t position = events[event_idx].first;
    if (position > previous_position && coverage >= 3) {
      if (!in_run) {
        run_start = previous_position;
        in_run = true;
      }
    } else if (position > previous_position && in_run) {
      loci.push_back({chromosome, run_start + 1, previous_position, cell});
      in_run = false;
    }

    while (event_idx < events.size() && events[event_idx].first == position) {
      coverage += events[event_idx].second;
      ++event_idx;
    }
    previous_position = position;
  }

  if (in_run) {
    loci.push_back({chromosome, run_start + 1, previous_position, cell});
  }
}

extern "C" SEXP multiomeR_bpcells_amulet_loci(
    SEXP external_pointer, SEXP selected_sexp) {
  if (TYPEOF(selected_sexp) != LGLSXP) {
    Rf_error("`selected` must be logical");
  }

  FragmentLoaderABI *loader = get_loader(external_pointer);
  int n_cells = Rf_length(selected_sexp);
  int *selected = LOGICAL(selected_sexp);
  for (int idx = 0; idx < n_cells; ++idx) {
    if (selected[idx] == NA_LOGICAL) {
      Rf_error("`selected` cannot contain missing values");
    }
  }

  std::vector<Locus> loci;
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> fragments_by_cell(
      n_cells);

  try {
    loader->restart();
    uint64_t fragments_seen = 0;
    while (loader->nextChr()) {
      uint32_t chromosome = loader->currentChr();
      for (auto &fragments : fragments_by_cell) {
        fragments.clear();
      }

      while (loader->load()) {
        uint32_t capacity = loader->capacity();
        uint32_t *cell = loader->cellData();
        uint32_t *start = loader->startData();
        uint32_t *end = loader->endData();
        for (uint32_t idx = 0; idx < capacity; ++idx) {
          if (cell[idx] >= static_cast<uint32_t>(n_cells)) {
            throw std::runtime_error("BPCells fragment cell index is out of range");
          }
          if (selected[cell[idx]] == TRUE) {
            fragments_by_cell[cell[idx]].push_back(
                std::make_pair(start[idx], end[idx]));
          }
        }
        fragments_seen += capacity;
        if ((fragments_seen & 0x3ffff) == 0) {
          R_CheckUserInterrupt();
        }
      }

      for (uint32_t cell = 0; cell < static_cast<uint32_t>(n_cells); ++cell) {
        if (selected[cell] == TRUE) {
          append_coverage_loci(
              fragments_by_cell[cell], chromosome, cell, loci);
        }
      }
    }
  } catch (const std::exception &exception) {
    Rf_error("%s", exception.what());
  }

  // scDblFinder splits fragments by cell before calculating coverage. Match
  // that cell-major return order exactly, including on multi-chromosome input.
  std::sort(
      loci.begin(),
      loci.end(),
      [](const Locus &left, const Locus &right) {
        if (left.cell != right.cell) {
          return left.cell < right.cell;
        }
        if (left.chr != right.chr) {
          return left.chr < right.chr;
        }
        if (left.start != right.start) {
          return left.start < right.start;
        }
        return left.end < right.end;
      });

  SEXP chromosome = PROTECT(Rf_allocVector(INTSXP, loci.size()));
  SEXP start = PROTECT(Rf_allocVector(INTSXP, loci.size()));
  SEXP end = PROTECT(Rf_allocVector(INTSXP, loci.size()));
  SEXP cell = PROTECT(Rf_allocVector(INTSXP, loci.size()));

  for (R_xlen_t idx = 0; idx < static_cast<R_xlen_t>(loci.size()); ++idx) {
    INTEGER(chromosome)[idx] = static_cast<int>(loci[idx].chr + 1);
    INTEGER(start)[idx] = static_cast<int>(loci[idx].start);
    INTEGER(end)[idx] = static_cast<int>(loci[idx].end);
    INTEGER(cell)[idx] = static_cast<int>(loci[idx].cell + 1);
  }

  SEXP output = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(output, 0, chromosome);
  SET_VECTOR_ELT(output, 1, start);
  SET_VECTOR_ELT(output, 2, end);
  SET_VECTOR_ELT(output, 3, cell);

  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, Rf_mkChar("chr"));
  SET_STRING_ELT(names, 1, Rf_mkChar("start"));
  SET_STRING_ELT(names, 2, Rf_mkChar("end"));
  SET_STRING_ELT(names, 3, Rf_mkChar("cell"));
  Rf_setAttrib(output, R_NamesSymbol, names);

  UNPROTECT(6);
  return output;
}
