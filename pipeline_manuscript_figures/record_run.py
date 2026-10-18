#!/usr/bin/env python3
"""Record the targets one build run completed, and fail if any errored.

Run from the checkout that built the store, after the run. Keeps the targets completed since
the run's build started (from its first driver log) and writes their runtimes, completion times,
bytes and errors to a TSV, with the run's node, which the benchmarks use to match each target
with the worker samples of that node. The resource and revision benchmarks read these files,
because later runs overwrite the store's metadata. The record also lists every file target
outside the store, such as configuration tables and Cell Ranger outputs, with in_store false:
targets records such an input file's modification time as its completion time, so a run that
first tracks it records no build of it. A file target that returns no file has no completion
time, so it takes the time at which a driver log reports it completed.

A target that errors after an earlier successful build keeps its old completion time in
the metadata, so errors and worker crashes are also read from the driver logs; each must
have been rebuilt without error since, as by a later rerun that this record spans.

Usage: record_run.py <store_meta> <output_tsv> <driver_log>...  (the record's driver logs in order)
"""
import csv
import datetime
import os
import re
import sys
import zoneinfo

meta_file, output_file, *driver_logs = sys.argv[1:]

# The driver logs local times and its node; the metadata records UTC days since 1970.
logs = [open(driver_log).readlines() for driver_log in driver_logs]
started = next(line for line in logs[0] if " Building script " in line)[:19]
node = re.search(r" on (\S+)", next(line for line in logs[0] if " Building " in line)).group(1)
lines = [line for log in logs for line in log]


def local_days(stamp):
    return datetime.datetime.fromisoformat(stamp).replace(tzinfo=zoneinfo.ZoneInfo("Europe/Copenhagen")).timestamp() / 86400


events = [(line[:22], *line.split()[2:5:2]) for line in lines
          if len(line.split()) > 4 and line.split()[3] in ("target", "branch")]
logged_errors = {name for _, event, name in events if event == "errored"}
completed_days = {name: local_days(stamp) for stamp, event, name in events if event == "completed"}
# A task whose worker crashed more often than crew allows gets no "errored" line.
logged_errors |= set(re.findall(r"crew worker of task '([^']+)' crashed", "".join(lines)))
start_days = local_days(started)

store = os.path.realpath(os.path.dirname(os.path.dirname(meta_file))) + os.sep
with open(meta_file) as meta:
    rows = {row["name"]: row for row in csv.DictReader(meta, delimiter="|")}
for row in rows.values():
    row["in_store"] = row["format"] != "file" or all(
        os.path.realpath(path).startswith(store) for path in row["path"].split("*") if path)
    if row["time"] == "t0s" and row["name"] in completed_days:
        row["time"] = f"t{completed_days[row['name']]}s"
completed = [row for row in rows.values() if row["type"] in ("stem", "branch")
             and (float(row["time"][1:-1]) >= start_days or not row["in_store"])]

with open(output_file, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(["name", "parent", "seconds", "completed_utc", "bytes", "error", "node", "in_store"])
    for row in sorted(completed, key=lambda row: row["time"]):
        days = float(row["time"][1:-1])
        writer.writerow([row["name"], row["parent"], row["seconds"],
                         datetime.datetime.fromtimestamp(days * 86400, datetime.timezone.utc).isoformat(timespec="seconds"),
                         row["bytes"], row["error"], node, str(row["in_store"]).upper()])

errored = [row["name"] for row in completed if row["error"]]
errored += sorted(name for name in logged_errors - set(errored)
                  if rows.get(name, {}).get("error") or float(rows.get(name, {"time": "t0s"})["time"][1:-1]) < start_days)
built = sum(row["in_store"] for row in completed)
print(f"{built} targets completed since {started} on {node}, and {len(completed) - built} input files; {len(errored)} errored.")
if errored or not built:
    sys.exit("Errored: " + ", ".join(errored[:20]) if errored else "No targets completed.")
