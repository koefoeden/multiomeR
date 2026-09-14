#!/usr/bin/env python3
"""Submit and report RC validation from the ordinary Esrum host."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

REPO = "koefoeden/multiomeR"
CONTEXTS = {
    "demo": ("full-demo/human", "full-demo/mouse"),
    "private": ("local-validation/differential-analysis",),
}
RC = re.compile(r"v\d+\.\d+\.\d+-rc\.\d+\Z")
ACTIVE = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "REQUEUED", "RESIZING"}


def command(*args):
    return subprocess.check_output(args, text=True).strip()


def save(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def publish(record, state, description):
    # Never interpolate private config, paths, exceptions, or logs into statuses.
    for context in CONTEXTS[record["profile"]]:
        command("gh", "api", "--method", "POST", f"repos/{REPO}/statuses/{record['sha']}",
                "-f", f"state={state}", "-f", f"context={context}",
                "-f", f"description={description}")


def terminal(job):
    # The installed watcher runs on the ordinary host with Slurm authentication.
    rows = command("sacct", "-n", "-P", "-j", job, "--format=JobIDRaw,State")
    for row in rows.splitlines():
        fields = row.split("|")
        if fields[0] == job and fields[1]:
            return fields[1].split()[0].rstrip("+") not in ACTIVE
    return False  # Accounting may lag. Never assume an absent job has finished.


def reconcile(record):
    result = Path(record["run"]) / "result.json"
    if result.exists():
        data = json.loads(result.read_text())
        if data.get("sha") != record["sha"] or data.get("profile_sha", "") != record["profile_sha"]:
            return "error", "Validation provenance mismatch; details retained locally"
        if data.get("exit_code") == 0 and data.get("validated") is True:
            return "success", "Pipeline execution and output checks passed"
        return "failure", "Validation failed; details retained locally"
    if not record.get("job"):
        return "error", "Submission interrupted; manual reconciliation required"
    if terminal(record["job"]):
        return "error", "Job ended without complete validation results; details retained locally"
    return "pending", "Validation queued or running on local compute"


def archive(repository, revision, destination):
    destination.mkdir(parents=True)
    bundle = destination.parent / (destination.name + ".tar")
    subprocess.run(["git", "-C", str(repository), "archive", revision, "-o", str(bundle)], check=True)
    subprocess.run(["tar", "-xf", str(bundle), "-C", str(destination)], check=True)
    bundle.unlink()


def submit(root, mirror, sha, tag, profile, profile_repository, records):
    key = f"{sha}/{profile}"
    attempt = records.get(key, {}).get("attempt", 0) + 1
    run = root / "runs" / f"{sha}-{profile}-{attempt}"
    run.mkdir(parents=True)
    checkout = run / "checkout"
    command("git", "--git-dir", str(mirror), "worktree", "add", "--detach", str(checkout), sha)
    profile_sha = ""
    if profile == "private":
        if command("git", "-C", str(profile_repository), "status", "--porcelain"):
            raise RuntimeError("Commit the private profile changes before validation")
        profile_sha = command("git", "-C", str(profile_repository), "rev-parse", "HEAD")
        archive(profile_repository, profile_sha, run / "profile")
    # Snapshot the harness so edits cannot change an already queued validation.
    shutil.copytree(Path(__file__).resolve().parent, run / "harness", ignore=shutil.ignore_patterns("__pycache__"))
    record = dict(sha=sha, tag=tag, profile=profile, profile_sha=profile_sha,
                  attempt=attempt, run=str(run), state="submitting")
    records[key] = record
    save(root / "state.json", records)
    publish(record, "pending", "Validation queued or running on local compute")
    job = command("sbatch", "--parsable", "--cpus-per-task=16", "--mem=256G",
                  "--time=24:00:00", "--export=NONE", "--job-name=multiomer-validation",
                  f"--output={run}/slurm.log", str(run / "harness/run.sh"),
                  str(checkout), str(root / "cache"), str(run), sha, profile_sha)
    record.update(job=job.split(";")[0], state="pending")
    save(root / "state.json", records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.home() / ".local/state/multiomer-validation")
    parser.add_argument("--profile-repository", type=Path,
                        default=os.environ.get("MULTIOMER_VALIDATION_PROFILE_REPO"))
    parser.add_argument("--retry", help="Retry SHA/demo or SHA/private after a terminal job")
    args = parser.parse_args()
    if not args.profile_repository:
        parser.error("Set --profile-repository or MULTIOMER_VALIDATION_PROFILE_REPO")
    os.umask(0o077)
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / "watch.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        mirror = root / "repository.git"
        if not mirror.exists():
            command("git", "clone", "--bare", f"https://github.com/{REPO}.git", str(mirror))
        command("git", "--git-dir", str(mirror), "fetch", "origin",
                "+refs/heads/main:refs/heads/main", "refs/tags/*:refs/tags/*")
        state_file = root / "state.json"
        records = json.loads(state_file.read_text()) if state_file.exists() else {}
        for record in records.values():
            if record["state"] in {"pending", "submitting"}:
                state, description = reconcile(record)
                if state != "pending" or record.get("description") != description:
                    publish(record, state, description)
                record.update(state=state, description=description)
                save(state_file, records)
        if args.retry:
            record = records[args.retry]
            if record["state"] not in {"success", "failure", "error"} or not record.get("job"):
                raise RuntimeError("Retry requires a reconciled terminal Slurm job")
            if not terminal(record["job"]):
                raise RuntimeError("The previous job is not confirmed terminal")
            submit(root, mirror, record["sha"], record["tag"], record["profile"], args.profile_repository, records)
            return
        tags = command("git", "--git-dir", str(mirror), "tag", "--list", "v*-rc.*").splitlines()
        for tag in tags:
            if not RC.fullmatch(tag):
                continue
            sha = command("git", "--git-dir", str(mirror), "rev-parse", f"{tag}^{{commit}}")
            command("git", "--git-dir", str(mirror), "merge-base", "--is-ancestor", sha, "main")
            if any(r["tag"] == tag and r["sha"] != sha for r in records.values()):
                raise RuntimeError(f"Release-candidate tag moved: {tag}")
            for profile in CONTEXTS:
                if f"{sha}/{profile}" not in records:
                    submit(root, mirror, sha, tag, profile, args.profile_repository, records)


if __name__ == "__main__":
    main()
