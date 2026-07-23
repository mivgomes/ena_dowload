#!/usr/bin/env python3
"""
Cross-check read counts for an ENA study against the submitted BAMs and the
ENA-generated FASTQs.

These runs were submitted as BAM. Each BAM holds collapsed (single-end) reads
plus some uncollapsed pairs, and ENA split that into three FASTQs per run:
  <run>.fastq.gz    collapsed reads
  <run>_1.fastq.gz  read 1 of the uncollapsed pairs
  <run>_2.fastq.gz  read 2 of the same pairs

ENA's read_count counts SPOTS: a collapsed read is 1 spot, a PAIR is also 1
spot. So the three sources should agree like this:

  read_count == collapsed + pairs        (both in the FASTQs and in the BAM)
  fastq <run>.fastq.gz  ==  BAM collapsed (unpaired) records
  fastq <run>_1.fastq.gz == fastq <run>_2.fastq.gz == BAM pairs

BAMs are named after the library, not the run, so we map run -> library via the
library_name column (or the submitted_ftp column when the table has one).

Usage:
  verify_read_counts.py --table input/PRJEB18629.tsv \
      --bam-dir input/submmitedBAMfiles \
      --fastq-dir input/generatedFASTQfiles \
      --out read_count_check.tsv

  # only a few runs, e.g. while debugging
  verify_read_counts.py ... --runs ERR1883466,ERR1883467

Notes:
  - Uses samtools and zcat when they are on PATH (much faster). Falls back to
    pure Python if not, so it still runs anywhere.
  - Secondary and supplementary BAM records are excluded from every count.
"""
import argparse
import csv
import gzip
import os
import shutil
import struct
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HAVE_SAMTOOLS = shutil.which("samtools") is not None
HAVE_ZCAT = shutil.which("zcat") is not None

def count_fastq(path):
    """Reads in a fastq(.gz) = lines / 4. Returns None if the file isn't there."""
    if not os.path.exists(path):
        return None
    if HAVE_ZCAT and path.endswith(".gz"):
        # zcat | wc -l is a lot faster than decompressing in Python
        zcat = subprocess.Popen(["zcat", path], stdout=subprocess.PIPE)
        wc = subprocess.run(["wc", "-l"], stdin=zcat.stdout, capture_output=True, text=True)
        zcat.stdout.close()
        if zcat.wait() != 0:
            raise RuntimeError(f"zcat failed on {path} (truncated or corrupt gzip?)")
        return int(wc.stdout.split()[0]) // 4
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n // 4


def _samtools_count(path, keep_flag, drop_flag):
    out = subprocess.run(
        ["samtools", "view", "-c", "-f", keep_flag, "-F", drop_flag, path],
        capture_output=True, text=True,
    )
    if out.returncode != 0:
        raise RuntimeError(f"samtools failed on {path}: {out.stderr.strip()}")
    return int(out.stdout.strip())


def _python_count(path):
    """Walk the BAM records ourselves and tally flags. Slow, but dependency-free."""
    collapsed = pairs = 0
    with gzip.open(path, "rb") as fh:
        if fh.read(4) != b"BAM\x01":
            raise RuntimeError(f"{path} is not a BAM file")
        l_text, = struct.unpack("<i", fh.read(4))
        fh.read(l_text)
        n_ref, = struct.unpack("<i", fh.read(4))
        for _ in range(n_ref):
            l_name, = struct.unpack("<i", fh.read(4))
            fh.read(l_name + 4)
        while True:
            head = fh.read(4)
            if len(head) < 4:
                break
            size, = struct.unpack("<i", head)
            block = fh.read(size)
            flag, = struct.unpack("<H", block[14:16])
            if flag & 0x900:            # secondary / supplementary
                continue
            if not flag & 0x1:          # not paired -> a collapsed read
                collapsed += 1
            elif flag & 0x40:           # first in pair -> count the pair once
                pairs += 1
    return collapsed, pairs


def count_bam(path):
    """(collapsed, pairs) for one BAM, ignoring secondary/supplementary records."""
    if not os.path.exists(path):
        return None, None
    if HAVE_SAMTOOLS:
        # 0x901 = paired | secondary | supplementary -> leaves primary unpaired reads
        collapsed = _samtools_count(path, "0", "0x901")
        pairs = _samtools_count(path, "0x40", "0x900")
        return collapsed, pairs
    return _python_count(path)


def bam_name_for(row):
    """BAMs are named after the library. Prefer submitted_ftp if the table has it."""
    ftp = (row.get("submitted_ftp") or "").strip()
    for part in ftp.split(";"):
        if part.endswith(".bam"):
            return os.path.basename(part)
    lib = (row.get("library_name") or "").strip()
    return f"{lib}.bam" if lib else ""


def safe(fn, *args):
    """Run a counter, but turn an unreadable/corrupt file into a note, not a crash.

    A broken file is exactly what we are hunting for, so one bad run must not
    take down the whole check.
    """
    try:
        return fn(*args), None
    except Exception as exc:
        return None, str(exc)


def check_run(row, bam_dir, fastq_dir):
    """Count one run from both sources and work out whether everything lines up."""
    run = row["run_accession"].strip()
    lib = (row.get("library_name") or "").strip()
    ena = row.get("read_count", "").strip()
    ena_n = int(ena) if ena.isdigit() else None

    notes = []
    merged, err = safe(count_fastq, os.path.join(fastq_dir, f"{run}.fastq.gz"))
    notes.append(err)
    r1, err = safe(count_fastq, os.path.join(fastq_dir, f"{run}_1.fastq.gz"))
    notes.append(err)
    r2, err = safe(count_fastq, os.path.join(fastq_dir, f"{run}_2.fastq.gz"))
    notes.append(err)

    bam_file = bam_name_for(row)
    if bam_file:
        counts, err = safe(count_bam, os.path.join(bam_dir, bam_file))
        notes.append(err)
        collapsed, pairs = counts if counts else (None, None)
    else:
        collapsed, pairs = None, None

    notes = [n for n in notes if n]

    # Treat an absent paired file as "this run had no uncollapsed reads".
    fq_spots = None if merged is None else merged + (r1 or 0)
    bam_spots = None if collapsed is None else collapsed + (pairs or 0)

    problems = ["UNREADABLE"] if notes else []
    if not any(os.path.exists(os.path.join(fastq_dir, f"{run}{s}.fastq.gz"))
               for s in ("", "_1", "_2")):
        problems.append("FASTQ_MISSING")
    if not bam_file or not os.path.exists(os.path.join(bam_dir, bam_file)):
        problems.append("BAM_MISSING")
    if r1 is not None and r2 is not None and r1 != r2:
        problems.append("PAIRS_UNEVEN")
    if ena_n is not None and fq_spots is not None and ena_n != fq_spots:
        problems.append("ENA_VS_FASTQ")
    if ena_n is not None and bam_spots is not None and ena_n != bam_spots:
        problems.append("ENA_VS_BAM")
    if merged is not None and collapsed is not None and merged != collapsed:
        problems.append("MERGED_VS_COLLAPSED")
    if r1 is not None and pairs is not None and r1 != pairs:
        problems.append("PAIRS_VS_BAM")

    return {
        "run_accession": run,
        "library_name": lib,
        "bam_file": bam_file,
        "ena_read_count": ena,
        "fq_merged": merged,
        "fq_r1": r1,
        "fq_r2": r2,
        "fq_spots": fq_spots,
        "bam_collapsed": collapsed,
        "bam_pairs": pairs,
        "bam_spots": bam_spots,
        "status": "OK" if not problems else ",".join(problems),
        "error": "; ".join(notes),
        "sample_description": (row.get("sample_description") or "").strip(),
    }


def main():
    ap = argparse.ArgumentParser(description="Compare ENA read_count with BAM and FASTQ counts.")
    ap.add_argument("--table", required=True, help="ENA report TSV (PRJEB18629.tsv)")
    ap.add_argument("--bam-dir", required=True, help="folder with the submitted BAMs")
    ap.add_argument("--fastq-dir", required=True, help="folder with the ENA-generated FASTQs")
    ap.add_argument("--out", default="read_count_check.tsv", help="where to write the full table")
    ap.add_argument("--runs", help="comma-separated run accessions (default: all rows)")
    ap.add_argument("--jobs", type=int, default=4, help="runs to count in parallel")
    args = ap.parse_args()

    with open(args.table) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("run_accession")]
    if args.runs:
        wanted = {r.strip() for r in args.runs.split(",")}
        rows = [r for r in rows if r["run_accession"].strip() in wanted]
    if not rows:
        sys.exit("No runs to check -- is --table the right file?")

    if not HAVE_SAMTOOLS:
        print("note: samtools not found, falling back to the Python BAM reader "
              "(slower on big files)\n", file=sys.stderr)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda r: check_run(r, args.bam_dir, args.fastq_dir), rows))
    results.sort(key=lambda d: d["run_accession"])

    def show(v):
        return "-" if v is None else str(v)

    header = f"{'run':13s} {'library':10s} {'ena':>10s} {'fq_merged':>10s} {'fq_r1':>7s} " \
             f"{'fq_spots':>10s} {'bam_coll':>10s} {'bam_pr':>7s} {'bam_spots':>10s} status"
    print(header)
    print("-" * len(header))
    for d in results:
        print(f"{d['run_accession']:13s} {d['library_name']:10s} {d['ena_read_count']:>10s} "
              f"{show(d['fq_merged']):>10s} {show(d['fq_r1']):>7s} {show(d['fq_spots']):>10s} "
              f"{show(d['bam_collapsed']):>10s} {show(d['bam_pairs']):>7s} "
              f"{show(d['bam_spots']):>10s} {d['status']}")

    cols = list(results[0].keys())
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(results)

    n_ok = sum(d["status"] == "OK" for d in results)
    print("-" * len(header))
    print(f"OK: {n_ok}   PROBLEM: {len(results) - n_ok}   (full table written to {args.out})")
    sys.exit(0 if n_ok == len(results) else 1)


if __name__ == "__main__":
    main()