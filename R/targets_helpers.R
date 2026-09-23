# Add packages to the default target packages.
w_def <- function(packages) {
  targets::tar_option_get("packages") %>%
    c(packages) %>%
    unique()
}
