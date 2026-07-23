# ena_dowload

## 1. /ena_download/download_from_ena.sh

Download the FASTQ files for an ENA study from the ENA filereport, with resume (-c) and md5 verification against ENA's published checksums.

Requires: curl, aria2, md5sum

```bash
conda create -n dl -c conda-forge aria2 -y
conda activate dl
```

For PRJEB18629 each run has three files on ENA:
   ERR*.fastq.gz     <- the generated merged file
   ERR*_1.fastq.gz   <- paired forward
   ERR*_2.fastq.gz   <- paired reverse

Default is to downlaod all, but you can use WHICH=merged or WHICH=paired if only interested on those

Usage:
 ```bash 
download_ena.sh
download_ena.sh PRJEB18629 input merged
```

## 2. /ena_download/verify_ena_downloads.py

Verify locally downloaded FASTQ files against ENA.

For each run it checks TWO independent things:
  1. md5   -> file integrity: local md5sum vs ENA's published fastq_md5
  2. reads -> local read count (gzip lines / 4) vs ENA's read_count field

A file that FAILS md5 is truncated/corrupt (re-download it) -- or was re-compressed after download (reads may be fine, but bytes differ from ENA).

Usage:
```bash
  verify_ena_downloads.py --accession PRJEB18629 --input-dir /input
  verify_ena_downloads.py --accession ERR1883403,ERR1883404 --input-dir /input
  verify_ena_downloads.py --accession PRJEB18629 --input-dir /input --skip-readcount
```

Notes:
  - --skip-readcount avoids decompressing every file (fast; md5 only).
  - ENA read_count counts SPOTS, not files. These runs were submitted as BAM and
    hold both collapsed and uncollapsed reads, so ENA emits three fastqs:
    <run>.fastq.gz (collapsed, 1 spot each) and <run>_1/_2.fastq.gz (a pair,
    also 1 spot). So the run reconciles as:
        read_count == reads(<run>.fastq.gz) + reads(<run>_1.fastq.gz)
    Comparing any single file against read_count always looks like a mismatch.



## 3. /ena_download/verify_read_counts.py

Cross-check read counts for an ENA study against the **submitted BAMs** and the
**ENA-generated FASTQs**.

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
```bash
# Assuming the dl enviroment was already created
# conda create -n dl -c conda-forge aria2 -y
# conda activate dl
conda install samtools

verify_read_counts.py --table input/PRJEB18629.tsv \
    --bam-dir input/submmitedBAMfiles \
    --fastq-dir input/generatedFASTQfiles \
    --out read_count_check.tsv
```

only a few runs, e.g. while debugging verify_read_counts.py ... --runs ERR1883466,ERR1883467

Notes:
- Uses samtools and zcat when they are on PATH (much faster). Falls back to pure Python if not, so it still runs anywhere.
- Secondary and supplementary BAM records are excluded from every count.