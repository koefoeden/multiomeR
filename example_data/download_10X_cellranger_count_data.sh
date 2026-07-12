#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="${1:-${script_dir}/public_core_cellranger_count_manifest.tsv}"
dest_root="${TENX_PUBLIC_DATA_DIR:-${script_dir}}"
overwrite="${OVERWRITE:-0}"

if ! command -v curl >/dev/null 2>&1; then
  echo "curl is required but was not found on PATH." >&2
  exit 1
fi

if [[ ! -f "${manifest}" ]]; then
  echo "Manifest not found: ${manifest}" >&2
  exit 1
fi

tail -n +2 "${manifest}" | while IFS=$'\t' read -r GEM_well_ID file_name url; do
  [[ -z "${GEM_well_ID}" ]] && continue

  output_dir="${dest_root}/${GEM_well_ID}/outs"
  output_file="${output_dir}/${file_name}"
  mkdir -p "${output_dir}"

  if [[ -s "${output_file}" && "${overwrite}" != "1" ]]; then
    echo "Skipping existing file: ${output_file}"
    continue
  fi

  tmp_file="$(mktemp "${output_file}.tmp.XXXXXX")"
  echo "Downloading ${url}"
  if curl --fail --location --retry 3 --retry-delay 5 --output "${tmp_file}" "${url}"; then
    mv "${tmp_file}" "${output_file}"
  else
    rm -f "${tmp_file}"
    exit 1
  fi
done
