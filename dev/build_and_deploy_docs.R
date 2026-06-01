unlink("docs", recursive = TRUE, force = TRUE)
status <- system2("quarto", c("render", "website"))
if (!identical(status, 0L)) {
  stop("Quarto render failed.", call. = FALSE)
}

system("bash dev/deploy_docs_to_gh_pages.sh")
