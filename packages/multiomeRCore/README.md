# multiomeRCore

`multiomeRCore` packages the generally reusable helper code that the multiomeR
pipeline sources directly from this directory. Ordinary multiomeR users do not
need to install this package: they can inspect and edit the files under `R/` as
normal pipeline source code.

Related standalone pipelines can instead install this subdirectory from an
immutable multiomeR Git commit. This provides versioned reuse without a second
copy of the helper implementations. Until the multiomeR repository is
published, installation from GitHub requires repository access and an
authenticated token.
