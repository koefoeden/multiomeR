group_split_by <- function(input, split_var) {
  grouped <- input %>%
    dplyr::group_by(.data[[split_var]])

  names <- dplyr::group_keys(grouped)[[split_var]] %>% as.character()

  split_w_names <- grouped %>%
    dplyr::group_split() %>%
    purrr::set_names(names) %>%
    sort_list()

  return(split_w_names)
}

# Input object works with character vectors and numeric vectors as well as factors
get_mixsorted_factor <- function(input_object) {
  # essentially re-convert to factor using a mixedsort order.
  reordered_factor <- input_object %>%
    as.character() %>%
    factor(., levels = gtools::mixedsort(unique(.)))
  return(reordered_factor)
}

sort_list <- function(list_input) {
  return(list_input[gtools::mixedorder(names(list_input))])
}

is_binary_vec <- function(vec) {
  all(vec %in% c(0, 1, NA))
}

is_non_binary_numeric_vec <- function(vec) {
  is.numeric(vec) && !is_binary_vec(vec)
}

is_count_matrix <- function(matrix) {
  sample_rows <- sample(nrow(matrix), min(10, nrow(matrix)))
  sample_cols <- sample(ncol(matrix), min(10, ncol(matrix)))
  all(matrix[sample_rows, sample_cols] %% 1 == 0)
}
