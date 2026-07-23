#!/usr/bin/env python3
"""
Verify locally downloaded FASTQ files against ENA.

For each run it checks TWO independent things:
  1. md5   -> file integrity: local md5sum vs ENA's published fastq_md5
  2. reads -> local read count (gzip lines / 4) vs ENA's read_count field

A file that FAILS md5 is truncated/corrupt (re-download it) -- or was re-compressed after download
(reads may be fine, but bytes differ from ENA).

Usage:
  verify_ena_downloads.py --accession PRJEB18629 --input-dir /input
  verify_ena_downloads.py --accession ERR1883403,ERR1883404 --input-dir /input
  verify_ena_downloads.py --accession PRJEB18629 --input-dir /input --skip-readcount

Notes:
  - --skip-readcount avoids decompressing every file (fast; md5 only).
  - ENA read_count counts SPOTS, not files. These runs were submitted as BAM and
    hold both collapsed and uncollapsed reads, so ENA emits three fastqs:
    <run>.fastq.gz (collapsed, 1 spot each) and <run>_1/_2.fastq.gz (a pair,
    also 1 spot). So the run reconciles as:
        read_count == reads(<run>.fastq.gz) + reads(<run>_1.fastq.gz)
    Comparing any single file against read_count always looks like a mismatch.
"""

import argparse
import gzip
import hashlib
import os
import sys
import urllib.request

ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
       "?accession={acc}&result=read_run"
       "&fields=run_accession,library_name,read_count,fastq_bytes,fastq_md5,fastq_ftp"
       "&format=tsv")

def fetch_report(accession):
    """Return list of dict rows from the ENA filereport for one accession."""
    url = ENA.format(acc=accession)
    with urllib.request.urlopen(url, timeout=60) as r:
        text = r.read().decode()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        sys.exit(f"ENA returned no records for '{accession}'. Check the accession.")
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]

def md5_of(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def count_reads(path):
    """Number of reads in a fastq(.gz): lines / 4."""
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n // 4


def main():
    ap = argparse.ArgumentParser(description="Verify FASTQ downloads against ENA.")
    ap.add_argument("--accession", required=True,
                    help="ENA study (PRJEB...) or comma-separated run accessions")
    ap.add_argument("--input-dir", required=True, help="folder holding the downloaded fastq(.gz)")
    ap.add_argument("--skip-readcount", action="store_true",
                    help="only check md5 (don't decompress to count reads)")
    args = ap.parse_args()

    rows = []
    for acc in args.accession.split(","):
        rows.extend(fetch_report(acc.strip()))

    ok = bad = missing = 0
    print(f"{'run':12s} {'library':10s} {'md5':10s} "
          f"{'collapsed':>11s} {'pairs':>8s} {'spots':>11s} {'ena_reads':>11s} status")
    print("-" * 95)
    for row in rows:
        run = row.get("run_accession", "?")
        lib = row.get("library_name", "")
        md5s = row.get("fastq_md5", "").split(";")
        ftps = row.get("fastq_ftp", "").split(";")
        ena_rc = row.get("read_count", "").strip()

        # md5 every file of the run; count collapsed and _1 separately so we can
        # rebuild ENA's spot count.
        n_files = n_md5_ok = 0
        collapsed = pairs = 0
        gone = False
        for md5_exp, ftp in zip(md5s, ftps):
            if not ftp:
                continue
            fname = os.path.basename(ftp)
            local = os.path.join(args.input_dir, fname)
            if not os.path.exists(local):
                print(f"{run:12s} {lib:10s} {'-':10s} "
                      f"{'-':>11s} {'-':>8s} {'-':>11s} {ena_rc:>11s} MISSING {fname}")
                gone = True
                continue
            n_files += 1
            n_md5_ok += (md5_of(local) == md5_exp)
            if args.skip_readcount:
                continue
            if fname.endswith("_2.fastq.gz"):
                continue                              # mate of _1, same spots
            elif fname.endswith("_1.fastq.gz"):
                pairs = count_reads(local)
            else:
                collapsed = count_reads(local)
        if gone:
            missing += 1
            continue

        md5_ok = (n_md5_ok == n_files)
        if args.skip_readcount:
            spots, reads_note = "-", ""
            cols = ("-", "-")
        else:
            spots = collapsed + pairs
            cols = (str(collapsed), str(pairs))
            reads_note = "" if (not ena_rc or str(spots) == ena_rc) else " READCOUNT_MISMATCH"
            spots = str(spots)
        status = "OK" if md5_ok else "MD5_FAIL"
        status += reads_note
        if md5_ok and not reads_note:
            ok += 1
        else:
            bad += 1
        print(f"{run:12s} {lib:10s} {f'{n_md5_ok}/{n_files} ok' if md5_ok else 'FAIL':10s} "
              f"{cols[0]:>11s} {cols[1]:>8s} {spots:>11s} {ena_rc:>11s} {status}")

    print("-" * 95)
    print(f"OK: {ok}   PROBLEM: {bad}   MISSING: {missing}")
    sys.exit(1 if (bad or missing) else 0)


if __name__ == "__main__":
    main()