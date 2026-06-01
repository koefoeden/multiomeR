#!/usr/bin/env bash

set -euo pipefail

mem_limit="${1:-256M}"
per_writer_mb="${2:-4096}"
parallel_writers="${3:-4}"
time_limit="${4:-00:30:00}"

usage() {
    cat <<'EOF'
Run a small Slurm job that tries to OOM from file-backed cache.

Usage:
    bash/test_file_cache_oom.sh [MEM_LIMIT] [PER_WRITER_MB] [PARALLEL_WRITERS] [TIME_LIMIT]

Examples:
    bash/test_file_cache_oom.sh
    bash/test_file_cache_oom.sh 128M 512 4 00:05:00
    bash/test_file_cache_oom.sh 256M 4096 4 00:30:00

Notes:
    - This test uses normal buffered writes with `dd if=/dev/zero`.
    - Buffered writes create dirty page-cache pages, which count toward the
      cgroup memory limit and are more likely to trigger a memcg OOM than
      read-only file cache.
    - The script prints cgroup usage together with cache/file, rss/anon,
      dirty, and writeback from `memory.stat`.
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    usage
    exit 0
fi

printf 'Launching test job with --mem=%s, per_writer_mb=%s, parallel_writers=%s, time=%s\n' \
    "$mem_limit" "$per_writer_mb" "$parallel_writers" "$time_limit"

job_output_file="$(mktemp)"
trap 'rm -f "$job_output_file"' EXIT

set +e
srun \
    --mem="$mem_limit" \
    --cpus-per-task=1 \
    --time="$time_limit" \
    --job-name="file-cache-oom-test" \
    bash -s -- "$per_writer_mb" "$parallel_writers" <<'EOF' 2>&1 | tee "$job_output_file"
set -euo pipefail

per_writer_mb="$1"
parallel_writers="$2"

is_cgroup_v2() {
    [[ -f /sys/fs/cgroup/cgroup.controllers ]]
}

read_value() {
    local path="$1"
    if [[ -r "$path" ]]; then
        tr -d '\n' < "$path"
    else
        printf 'NA'
    fi
}

bytes_to_mib() {
    local bytes="$1"
    awk -v value="$bytes" 'BEGIN {
        if (value == "" || value == "NA") {
            printf "NA"
        } else {
            printf "%.1f", value / 1024 / 1024
        }
    }'
}

job_root_from_self() {
    local cg_rel
    if is_cgroup_v2; then
        cg_rel="$(awk -F: '$1 == "0" { print $3 }' /proc/self/cgroup)"
        printf '%s\n' "$cg_rel" | sed -E 's#/step_[^/]+/task_[^/]+$##; s#/task_[^/]+$##'
    else
        cg_rel="$(awk -F: '$2 == "memory" { print $3 }' /proc/self/cgroup)"
        printf '%s\n' "$cg_rel" | sed -E 's#/step_[^/]+/task_[^/]+$##; s#/task_[^/]+$##'
    fi
}

scope_path_from_self() {
    local cg_root
    local job_rel
    local job_path
    if is_cgroup_v2; then
        cg_root="/sys/fs/cgroup"
    else
        cg_root="/sys/fs/cgroup/memory"
    fi

    job_rel="$(job_root_from_self)"
    job_path="${cg_root}${job_rel}"

    if [[ -d "${job_path}/step_batch" ]]; then
        printf '%s\n' "${job_path}/step_batch"
    else
        printf '%s\n' "$job_path"
    fi
}

memory_stat_value() {
    local key="$1"
    local stat_file="$2"
    awk -v wanted="$key" '$1 == wanted { print $2; found = 1 } END { if (!found) print "NA" }' "$stat_file" 2>/dev/null
}

print_snapshot() {
    local label="$1"
    local scope_path="$2"
    local usage_file peak_file limit_file stat_file
    local usage peak limit cache rss dirty writeback

    if is_cgroup_v2; then
        usage_file="memory.current"
        peak_file="memory.peak"
        limit_file="memory.max"
        stat_file="${scope_path}/memory.stat"
        cache="$(memory_stat_value file "$stat_file")"
        rss="$(memory_stat_value anon "$stat_file")"
        dirty="$(memory_stat_value file_dirty "$stat_file")"
        writeback="$(memory_stat_value file_writeback "$stat_file")"
    else
        usage_file="memory.usage_in_bytes"
        peak_file="memory.max_usage_in_bytes"
        limit_file="memory.limit_in_bytes"
        stat_file="${scope_path}/memory.stat"
        cache="$(memory_stat_value total_cache "$stat_file")"
        rss="$(memory_stat_value total_rss "$stat_file")"
        dirty="$(memory_stat_value total_dirty "$stat_file")"
        writeback="$(memory_stat_value total_writeback "$stat_file")"
    fi

    usage="$(read_value "${scope_path}/${usage_file}")"
    peak="$(read_value "${scope_path}/${peak_file}")"
    limit="$(read_value "${scope_path}/${limit_file}")"

    printf '%s usage_mib=%s peak_mib=%s limit_mib=%s cache_mib=%s rss_mib=%s dirty_mib=%s writeback_mib=%s\n' \
        "$label" \
        "$(bytes_to_mib "$usage")" \
        "$(bytes_to_mib "$peak")" \
        "$(bytes_to_mib "$limit")" \
        "$(bytes_to_mib "$cache")" \
        "$(bytes_to_mib "$rss")" \
        "$(bytes_to_mib "$dirty")" \
        "$(bytes_to_mib "$writeback")"
}

scope_path="$(scope_path_from_self)"
work_root="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
work_dir="${work_root}/file-cache-oom-${SLURM_JOB_ID:-$$}"
mkdir -p "$work_dir"
trap 'rm -rf "$work_dir"' EXIT

printf 'host=%s job_id=%s scope_path=%s work_dir=%s\n' \
    "$(hostname)" \
    "${SLURM_JOB_ID:-NA}" \
    "$scope_path" \
    "$work_dir"
printf 'SLURM_TEST_JOB_ID=%s\n' "${SLURM_JOB_ID:-NA}"

print_snapshot "before" "$scope_path"
printf 'Starting %s buffered writers, each writing %s MiB\n' "$parallel_writers" "$per_writer_mb"

declare -a pids=()
for writer_idx in $(seq 1 "$parallel_writers"); do
    dd if=/dev/zero of="${work_dir}/writer_${writer_idx}.bin" bs=1M count="$per_writer_mb" status=none &
    pids+=("$!")
done

while true; do
    any_running=0
    for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            any_running=1
            break
        fi
    done

    print_snapshot "during" "$scope_path"

    if ((any_running == 0)); then
        break
    fi
    sleep 1
done

for pid in "${pids[@]}"; do
    wait "$pid"
done

print_snapshot "after_writes" "$scope_path"

printf 'Reading the files back through /dev/null to show clean file cache state\n'
for file_path in "${work_dir}"/writer_*.bin; do
    dd if="$file_path" of=/dev/null bs=16M status=none
done

print_snapshot "after_reads" "$scope_path"
printf 'Done\n'
EOF
srun_exit_code="${PIPESTATUS[0]}"
set -e

job_id="$(awk -F= '/^SLURM_TEST_JOB_ID=/{print $2; exit}' "$job_output_file")"

if [[ -n "${job_id:-}" ]] && [[ "$job_id" != "NA" ]]; then
    printf '\nSlurm accounting for job %s\n' "$job_id"
    sacct_output=""
    for _ in $(seq 1 10); do
        sacct_output="$(sacct -j "$job_id" --format=JobID,JobName%30,State,ExitCode,ReqMem,MaxRSS,MaxRSSNode,MaxRSSTask -P 2>/dev/null || true)"
        if [[ -n "$sacct_output" ]] && [[ "$(printf '%s\n' "$sacct_output" | awk 'END { print NR }')" -gt 1 ]]; then
            break
        fi
        sleep 2
    done

    if [[ -n "$sacct_output" ]]; then
        printf '%s\n' "$sacct_output"
    else
        printf 'sacct did not return accounting rows yet for job %s\n' "$job_id"
    fi
else
    printf '\nCould not determine the Slurm job id from srun output, so MaxRSS could not be queried.\n'
fi

exit "$srun_exit_code"
