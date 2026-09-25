# multiomeRCore

`multiomeRCore` packages the generally reusable helper code that the multiomeR
pipeline sources directly from this directory. Ordinary multiomeR users do not
need to install this package: they can inspect and edit the files under `R/` as
normal pipeline source code.

Related standalone pipelines can instead install this subdirectory from an
immutable multiomeR Git commit. This provides versioned reuse without a second
copy of the helper implementations.

`save_plots_structured()` renders replacements before publishing them. Each
target records its image paths under the store's `plot_inventory/`; later saves
remove superseded owned files, including for empty results. Keep this directory
with the store. Files belonging to sibling branches or unrelated tools are
preserved. For outputs created before inventories existed, cleanup can recover
paths registered in targets; it leaves unregistered historical files alone
rather than guess ownership.

Saved ggplot and patchwork captions are left-aligned. Configuration values added
with `add_plot_parameters()` are shortened to 40 characters by default.

Changes here affect consumers sourcing these files immediately. Consumers using
an installed, pinned package need an updated package revision to receive them.
