project_root <- rprojroot::find_root(
  rprojroot::has_file("pixi.toml"),
  path = getwd()
)
setwd(project_root)

testthat::test_dir(file.path(project_root, "tests", "testthat"))
