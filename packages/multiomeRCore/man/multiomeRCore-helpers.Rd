\name{multiomeRCore-helpers}
\alias{multiomeRCore-helpers}
\alias{apply_crew_controller_options}
\alias{get_structured_file_path}
\alias{get_tar_resources}
\alias{run_w_error_check}
\alias{save_plots_structured}
\alias{tar_name_wo_suffixes}
\title{Reusable helpers for targets pipelines}
\description{
General helpers for process execution, Crew resource selection, and structured
target output paths. multiomeR sources the underlying R files directly, while
standalone pipelines can install and attach this package.
}
\usage{
apply_crew_controller_options(controller_setup)

get_structured_file_path(filetype = NULL, override_suffix = NULL)

get_tar_resources(
  cores_req = NULL,
  RAM_GB_req = NULL,
  RAM_GB_per_extra_core = NULL,
  gpus_req = 0,
  controller_resources = NULL
)

run_w_error_check(command_string, arguments_chr = character(), ...)

save_plots_structured(
  plots,
  filetype = "png",
  override_suffix = NULL,
  dyn_suffix_in_subdir = FALSE,
  target_name = targets::tar_name(),
  discrete_legend_n_distinct_max = 20,
  save_serialized_plot_objects = multiomeR_save_serialized_plot_objects,
  ...
)

tar_name_wo_suffixes(target_name = targets::tar_name())
}
\arguments{
  \item{controller_setup}{A named list containing a controller resource table
  and corresponding Crew controllers.}
  \item{filetype}{File extension without a leading dot.}
  \item{override_suffix}{Optional replacement for the target-derived suffix.}
  \item{cores_req}{Minimum requested CPU cores.}
  \item{RAM_GB_req}{Minimum requested RAM in GB.}
  \item{RAM_GB_per_extra_core}{Additional RAM per requested core beyond the first.}
  \item{gpus_req}{Minimum requested GPUs.}
  \item{controller_resources}{Optional controller resource table.}
  \item{command_string}{Executable name or path.}
  \item{arguments_chr}{Character vector of command-line arguments.}
  \item{plots}{A plot object or non-empty list of plot objects.}
  \item{dyn_suffix_in_subdir}{Whether to place a dynamic branch suffix in a subdirectory.}
  \item{target_name}{Target name used to derive the output path.}
  \item{discrete_legend_n_distinct_max}{Maximum number of discrete legend entries to retain.}
  \item{save_serialized_plot_objects}{Whether to save mirrored serialized plot objects.}
  \item{...}{Additional arguments forwarded to the process or graphics helper.}
}
\value{
The configured resource object, command output path, structured path, saved
plot path, or target stem appropriate to the selected helper.
}
