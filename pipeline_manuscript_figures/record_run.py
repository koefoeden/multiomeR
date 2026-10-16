#!/usr/bin/env python3
"""Record the targets one build run completed, and fail if any errored.

Reads the store metadata after a run of pipeline_manuscript_docs/scripts/run_dev_local.sh
has copied the store back, keeps the targets completed since the run's build started (from
its driver log), and writes their runtimes and completion times to a TSV. The resource and
revision benchmarks read these files, because later runs overwrite the store's metadata.

Usage: record_run.py <driver_log> <store_meta> <output_tsv>
"""
import csv
import datetime
import sys
import zoneinfo

driver_log, meta_file, output_file = sys.argv[1:]

# run_dev_local.sh logs local times; the metadata records UTC days since 1970.
with open(driver_log) as log:
    started = next(line for line in log if " Building script " in line)[:19]
start_days = (datetime.datetime.fromisoformat(started).replace(tzinfo=zoneinfo.ZoneInfo("Europe/Copenhagen"))
              .timestamp() / 86400)

with open(meta_file) as meta:
    rows = {row["name"]: row for row in csv.DictReader(meta, delimiter="|")}
completed = [row for row in rows.values()
             if row["type"] in ("stem", "branch") and float(row["time"][1:-1]) >= start_days]

with open(output_file, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(["name", "parent", "seconds", "completed_utc", "bytes", "error"])
    for row in sorted(completed, key=lambda row: row["time"]):
        days = float(row["time"][1:-1])
        writer.writerow([row["name"], row["parent"], row["seconds"],
                         datetime.datetime.fromtimestamp(days * 86400, datetime.timezone.utc).isoformat(timespec="seconds"),
                         row["bytes"], row["error"]])

errored = [row["name"] for row in completed if row["error"]]
print(f"{len(completed)} targets completed since {started}; {len(errored)} errored.")
if errored or not completed:
    sys.exit("Errored: " + ", ".join(errored[:20]) if errored else "No targets completed.")
