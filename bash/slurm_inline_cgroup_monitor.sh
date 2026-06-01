#!/usr/bin/env bash

set -uo pipefail

interval_s="${1:-${SLURM_CGROUP_MONITOR_INTERVAL:-30}}"
tag="${SLURM_CGROUP_MONITOR_TAG:-CGROUP_MEM}"

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

bytes_to_gib() {
    local bytes="$1"

    awk -v value="$bytes" 'BEGIN {
        if (value == "" || value == "NA") {
            printf "NA"
        } else {
            printf "%.1f", value / 1024 / 1024 / 1024
        }
    }'
}

percent_used() {
    local used_bytes="$1"
    local limit_bytes="$2"

    awk -v used="$used_bytes" -v limit="$limit_bytes" 'BEGIN {
        if (used == "" || used == "NA" || limit == "" || limit == "NA" || limit <= 0 || limit > 9e18) {
            printf "NA"
        } else {
            printf "%.1f", 100 * used / limit
        }
    }'
}

status_from_metrics() {
    local used_pct="$1"
    local failcnt="$2"
    local oom_kill="$3"
    local limit_bytes="$4"

    if [[ "$oom_kill" != "NA" ]] && (( oom_kill > 0 )); then
        printf '%s\n' "OOM_HIT"
        return
    fi

    if [[ "$failcnt" != "NA" ]] && (( failcnt > 0 )); then
        printf '%s\n' "LIMIT_HIT"
        return
    fi

    awk -v pct="$used_pct" -v limit="$limit_bytes" 'BEGIN {
        if (limit == "" || limit == "NA" || limit > 9e18) {
            print "NO_LIMIT"
        } else if (pct == "" || pct == "NA") {
            print "UNKNOWN"
        } else if (pct >= 95) {
            print "DANGER"
        } else if (pct >= 85) {
            print "WARN"
        } else {
            print "OK"
        }
    }'
}

resolve_scope_path() {
    local cg_root
    local cg_rel
    local job_rel
    local job_path
    local scope_name
    local scope_path

    if is_cgroup_v2; then
        cg_root="/sys/fs/cgroup"
        cg_rel="$(awk -F: '$1 == "0" { print $3 }' /proc/self/cgroup)"
    else
        cg_root="/sys/fs/cgroup/memory"
        cg_rel="$(awk -F: '$2 == "memory" { print $3 }' /proc/self/cgroup)"
    fi

    job_rel="$(printf '%s\n' "$cg_rel" | sed -E 's#/step_[^/]+/task_[^/]+$##; s#/task_[^/]+$##')"
    job_path="${cg_root}${job_rel}"
    scope_name="job"
    scope_path="$job_path"

    if [[ -d "${job_path}/step_batch" ]]; then
        scope_name="batch"
        scope_path="${job_path}/step_batch"
    fi

    printf '%s|%s\n' "$scope_name" "$scope_path"
}

collect_metrics() {
    local scope_name="$1"
    local scope_path="$2"
    local usage_file
    local peak_file
    local limit_file
    local events_file
    local failcnt
    local under_oom
    local oom_kill

    if is_cgroup_v2; then
        usage_file="memory.current"
        peak_file="memory.peak"
        limit_file="memory.max"
        events_file="memory.events"
    else
        usage_file="memory.usage_in_bytes"
        peak_file="memory.max_usage_in_bytes"
        limit_file="memory.limit_in_bytes"
        events_file="memory.oom_control"
    fi

    usage="$(read_value "${scope_path}/${usage_file}")"
    peak="$(read_value "${scope_path}/${peak_file}")"
    limit="$(read_value "${scope_path}/${limit_file}")"

    if is_cgroup_v2; then
        failcnt="$(awk '$1 == "max" { print $2 }' "${scope_path}/${events_file}" 2>/dev/null || printf 'NA')"
        under_oom="$(awk '$1 == "oom" { print $2 }' "${scope_path}/${events_file}" 2>/dev/null || printf 'NA')"
        oom_kill="$(awk '$1 == "oom_kill" { print $2 }' "${scope_path}/${events_file}" 2>/dev/null || printf 'NA')"
    else
        failcnt="$(read_value "${scope_path}/memory.failcnt")"
        under_oom="$(awk '$1 == "under_oom" { print $2 }' "${scope_path}/${events_file}" 2>/dev/null || printf 'NA')"
        oom_kill="$(awk '$1 == "oom_kill" { print $2 }' "${scope_path}/${events_file}" 2>/dev/null || printf 'NA')"
    fi

    printf '%s|%s|%s|%s|%s|%s|%s\n' \
        "$scope_name" \
        "$usage" \
        "$peak" \
        "$limit" \
        "$failcnt" \
        "$under_oom" \
        "$oom_kill"
}

job_display_id="${SLURM_JOB_ID:-NA}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    job_display_id="${job_display_id}_${SLURM_ARRAY_TASK_ID}"
fi

printf '[%s] ts=%s host=%s job_id=%s interval_s=%s status=START\n' \
    "$tag" \
    "$(date '+%F %T')" \
    "$(hostname)" \
    "$job_display_id" \
    "$interval_s"

while true; do
    IFS='|' read -r scope_name scope_path <<< "$(resolve_scope_path)"
    IFS='|' read -r scope_name current_bytes peak_bytes limit_bytes failcnt under_oom oom_kill <<< "$(collect_metrics "$scope_name" "$scope_path")"

    current_gb="$(bytes_to_gib "$current_bytes")"
    peak_gb="$(bytes_to_gib "$peak_bytes")"
    limit_gb="$(bytes_to_gib "$limit_bytes")"
    used_pct="$(percent_used "$current_bytes" "$limit_bytes")"
    status="$(status_from_metrics "$used_pct" "$failcnt" "$oom_kill" "$limit_bytes")"

    printf '[%s] ts=%s host=%s job_id=%s scope=%s cur_gb=%s peak_gb=%s lim_gb=%s used_pct=%s failcnt=%s under_oom=%s oom_kill=%s status=%s\n' \
        "$tag" \
        "$(date '+%F %T')" \
        "$(hostname)" \
        "$job_display_id" \
        "$scope_name" \
        "$current_gb" \
        "$peak_gb" \
        "$limit_gb" \
        "$used_pct" \
        "$failcnt" \
        "$under_oom" \
        "$oom_kill" \
        "$status"

    sleep "$interval_s"
done
