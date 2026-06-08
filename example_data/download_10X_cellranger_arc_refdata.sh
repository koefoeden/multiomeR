#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="${1:-${script_dir}/public_core_cellranger_arc_refdata_manifest.tsv}"
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

mkdir -p "${dest_root}"

tail -n +2 "${manifest}" | while IFS=$'\t' read -r archive_name extract_dir url; do
  [[ -z "${archive_name}" ]] && continue

  archive_file="${dest_root}/${archive_name}"
  extract_path="${dest_root}/${extract_dir}"

  if [[ -s "${archive_file}" && "${overwrite}" != "1" ]]; then
    echo "Skipping existing archive: ${archive_file}"
  else
    tmp_file="$(mktemp "${archive_file}.tmp.XXXXXX")"
    echo "Downloading ${url}"
    if curl --fail --location --retry 3 --retry-delay 5 --output "${tmp_file}" "${url}"; then
      mv "${tmp_file}" "${archive_file}"
    else
      rm -f "${tmp_file}"
      exit 1
    fi
  fi

  if [[ -d "${extract_path}" && "${overwrite}" != "1" ]]; then
    echo "Skipping existing reference directory: ${extract_path}"
  else
    echo "Extracting ${archive_file}"
    tar -xzf "${archive_file}" -C "${dest_root}"
  fi
done
