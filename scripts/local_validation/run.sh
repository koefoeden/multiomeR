#!/usr/bin/env bash
# Run only in an isolated checkout on an allocated compute node.
set -euo pipefail
umask 077
checkout=$(realpath "${1:?checkout}")
cache=$(realpath -m "${2:?cache directory}")
run=$(realpath "${3:?run directory}")
sha=${4:?commit SHA}
profile_sha=${5:-}
harness="$run/harness"
export PATH="$HOME/.pixi/bin:$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin:$PATH"
# The standard heavy Crew workers use six cores. Do not inherit the launcher's
# OpenMP cap (which may describe a smaller interactive allocation).
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=1
mkdir -p "$cache" "$run"
finish() {
  code=$?
  trap - EXIT
  validated=false
  [[ $code == 0 && -f "$run/validated" ]] && validated=true
  printf '{"sha":"%s","profile_sha":"%s","exit_code":%d,"validated":%s}\n' "$sha" "$profile_sha" "$code" "$validated" > "$run/result.json.tmp"
  mv "$run/result.json.tmp" "$run/result.json"
  exit "$code"
}
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
cd "$checkout"
test "$(git rev-parse HEAD)" = "$sha"
test -z "$(git status --porcelain)"
test ! -e outputs/meta/meta
test ! -e "$run/validated"
if [[ -n "$profile_sha" ]]; then
  test -f "$run/profile/validate.R"
  # Fixed destinations only; private files never enter the public source checkout.
  cp "$run/profile/cfg_GEM_wells.tsv" cfg_GEM_wells.tsv
  cp "$run/profile/cfg_aggregations.yaml" cfg_aggregations.yaml
  cp "$run/profile/module_differential_analyses/cfg.yaml" module_differential_analyses/cfg.yaml
  if [[ -d "$run/profile/example_data" ]]; then
    cp -R "$run/profile/example_data/." example_data/
  fi
  printf "{}\n" > module_genetic_enrichment/cfg.yaml
fi
env_key=$(cat pixi.toml pixi.lock scripts/install_r_github_packages.R | sha256sum | cut -d' ' -f1)
mkdir -p "$cache/environments/$env_key"
ln -s "$cache/environments/$env_key" .pixi
(
  flock 9
  pixi run --locked --use-environment-activation-cache install-r-github-packages
) 9> "$cache/environments/$env_key.lock"
mkdir -p "$cache/AnnotationHub" outputs/files
ln -s "$cache/AnnotationHub" outputs/files/AnnotationHub
if [[ -z "$profile_sha" ]]; then
  pixi run --locked --use-environment-activation-cache Rscript "$harness/prepare.R" "$cache" "$run"
  (
    flock 9
    TENX_PUBLIC_DATA_DIR="$(cat "$run/input-root")" bash example_data/download_10X_cellranger_count_data.sh "$run/manifest.tsv"
  ) 9> "$cache/inputs.lock"
fi
export MULTIOMER_VALIDATION_RUN="$run"
pixi run --locked --use-environment-activation-cache Rscript "$harness/validate.R" &
child=$!
trap 'kill -TERM "$child" 2>/dev/null || true; wait "$child" || true; exit 143' TERM
trap 'kill -INT "$child" 2>/dev/null || true; wait "$child" || true; exit 130' INT
wait "$child"
