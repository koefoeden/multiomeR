# Version 0.5.0 preparation validation

Checked on 2026-09-11 in an isolated worktree, using the configured Pixi
runtime with `--frozen --no-install --use-environment-activation-cache`. The
installed environment was reused without installing or updating packages.
This does not establish a fresh-install test.

- Default public manifest: 428 targets. Explicit peak-gene module opt-in:
  463 targets, including 35 peak-gene targets. No legacy annotation method or
  labeller remains in target commands. The default demo configuration is retained.
- Focused QC, gene-feature, cluster-marker and GWAS-layout checks pass.
- Numerical scoring and peak-gene tests: 91 assertions, zero failures, warnings
  or skips. This includes signed observed/control UCell means and marker-deletion
  variants against UCell 2.14.0, chunk/worker agreement, HC3 parity and compact
  global-FDR reconstruction. The R process exited successfully; R printed a
  parallel-child cleanup message at shutdown, after the successful test result.
- A synthetic 2,000-gene, 600-cell signed-marker check with 999 controls and two
  workers completed its summary in 0.99 seconds; the summary occupied 4.55 MiB.
  The entire smoke-check process, including startup, evidence plotting and
  SuSiE, peaked at 466 MiB RSS. These small-fixture observations do not predict
  production runtime or memory for larger marker panels and cluster counts.
- A nonempty four-peak, 120-aggregate SuSiE smoke fit converged and assigned PIP 1
  to its planted strong signal. This checks execution and a basic expected
  numerical behavior, not causal inference or biological calibration.
- Manifest-derived diagrams, both Quarto books and the LLM Markdown export
  were regenerated. Whitespace and new-public-material checks passed.

Checks use synthetic data and in-process manifests. No production targets store
was read, copied, changed or run. No full public-demo or cohort pipeline replay
was performed. Real-data annotation review and fresh-install CI remain distinct
release acceptance evidence; do not describe these checks as biological validation.

Main commands (from this worktree):

```sh
pixi run --frozen --no-install --use-environment-activation-cache Rscript -e 'testthat::test_dir("tests/testthat", filter="scoring-parity|peak-gene-sidecars", stop_on_failure=TRUE)'
pixi run --frozen --no-install --use-environment-activation-cache Rscript website/figures/human_curated/graphs_v2.R
pixi run --frozen --no-install --use-environment-activation-cache quarto render website
pixi run --frozen --no-install --use-environment-activation-cache quarto render website/implementation
pixi run --frozen --no-install --use-environment-activation-cache -e dev export-website-llm-markdown
git diff --check
```
