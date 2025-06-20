#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path
import pandas as pd
from collections import defaultdict
from shutil import which
from datetime import datetime

def run_command(command, log_path=None, error_message="Command failed"):
    print(f"\nRunning: {command}")
    with subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
        with open(log_path, "a") if log_path else open(os.devnull, "w") as log_file:
            for line in proc.stdout:
                print(line, end="")
                log_file.write(line)
        if proc.wait() != 0:
            raise RuntimeError(error_message)

def find_executable(tool_name, fallback_path=None):
    if which(tool_name):
        return tool_name
    elif fallback_path and Path(fallback_path).exists():
        return str(fallback_path)
    else:
        raise FileNotFoundError(f"{tool_name} not found in PATH or provided path.")

def main():
    parser = argparse.ArgumentParser(description="ARIMA Genomics Mapping Pipeline with Tech/Bio Replicate Handling")
    parser.add_argument("fastq_map", help="TSV file: read1 read2 tech_id bio_id sample_label")
    parser.add_argument("prefix", help="Output file prefix")
    parser.add_argument("reference", help="Reference FASTA (indexed)")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("threads", type=int, help="Number of threads")
    parser.add_argument("--filter", default="filter_five_end.pl", help="Path to 5' filter script")
    parser.add_argument("--combiner", default="two_read_bam_combiner.pl", help="Path to two-read combiner script")
    parser.add_argument("--stats", default="bam_stats.pl", help="Path to Perl BAM stats script")
    args = parser.parse_args()

    filter_script = find_executable("filter_five_end.pl", args.filter)
    combiner_script = find_executable("two_read_bam_combiner.pl", args.combiner)
    stats_script = find_executable("bam_stats.pl", args.stats)

    bwa = "bwa"
    samtools = "samtools"
    picard = "picard"

    # Setup output dirs
    outdir = Path(args.outdir)
    dirs = {
        "raw": outdir / "raw",
        "filtered": outdir / "filtered",
        "combined": outdir / "combined",
        "paired": outdir / "paired",
        "replicate": outdir / "replicate",
        "merged": outdir / "merged",
        "temp": outdir / "temp",
        "logs": outdir / "logs"
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.fastq_map, sep="\t", header=None,
                     names=["read1", "read2", "tech_id", "bio_id", "label"])
    sample_label = df["label"].unique()[0]
    assert len(df["label"].unique()) == 1, "Multiple sample labels found; expected only one."

    tech_groups = defaultdict(list)
    bio_groups = defaultdict(list)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = dirs["logs"] / f"{args.prefix}_{timestamp}.log"

    for _, row in df.iterrows():
        r1, r2, tech_id, bio_id, label = row
        raw1 = dirs["raw"] / f"{tech_id}_1.bam"
        raw2 = dirs["raw"] / f"{tech_id}_2.bam"
        filt1 = dirs["filtered"] / f"{tech_id}_1.bam"
        filt2 = dirs["filtered"] / f"{tech_id}_2.bam"
        combined = dirs["combined"] / f"{tech_id}.bam"
        rg_bam = dirs["paired"] / f"{tech_id}.bam"

        print(f"\n### Mapping: {tech_id} ({bio_id})")
        run_command(f"{bwa} mem -t{args.threads} {args.reference} {r1} | {samtools} view -Sb - > {raw1}", log_file)
        run_command(f"{samtools} view -h {raw1} | perl {filter_script} | {samtools} view -@{args.threads} -Sb - > {filt1}", log_file)
        raw1.unlink()

        run_command(f"{bwa} mem -t{args.threads} {args.reference} {r2} | {samtools} view -Sb - > {raw2}", log_file)
        run_command(f"{samtools} view -h {raw2} | perl {filter_script} | {samtools} view -@{args.threads} -Sb - > {filt2}", log_file)
        raw2.unlink()

        run_command(f"perl {combiner_script} {filt1} {filt2} | {samtools} view -@{args.threads} -Sb > {combined}", log_file)
        run_command(f"{picard} AddOrReplaceReadGroups I={combined} O={rg_bam} ID={tech_id} LB={tech_id} SM={label} PL=ILLUMINA PU=none", log_file)

        tech_groups[tech_id].append(str(rg_bam))
        bio_groups[bio_id].append(str(rg_bam))

    for tech_id, bams in tech_groups.items():
        print(f"\n### Technical Replicate: {tech_id}")
        if len(bams) > 1:
            merged = dirs["replicate"] / f"{tech_id}.bam"
            inputs = ' '.join([f"I={bam}" for bam in bams])
            run_command(f"{picard} MergeSamFiles {inputs} O={merged} USE_THREADING=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT", log_file)
        else:
            merged = Path(bams[0])

        dedup = dirs["replicate"] / f"{tech_id}.dedup.bam"
        run_command(
            f"{picard} MarkDuplicates I={merged} O={dedup} M={dirs['replicate']}/metrics.{tech_id}.txt "
            f"TMP_DIR={dirs['temp']} ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true", log_file
        )
        run_command(f"{samtools} index {dedup}", log_file)
        run_command(f"perl {stats_script} {dedup} > {dirs['replicate']}/{tech_id}.dedup.bam.stats", log_file)

        for bio_id in bio_groups:
            bio_groups[bio_id] = [str(dedup) if str(b).endswith(f"{tech_id}.bam") else b for b in bio_groups[bio_id]]

    for bio_id, bams in bio_groups.items():
        print(f"\n### Biological Replicate: {bio_id}")
        if len(bams) > 1:
            bio_merged = dirs["merged"] / f"{bio_id}.bam"
            inputs = ' '.join([f"I={bam}" for bam in bams])
            run_command(f"{picard} MergeSamFiles {inputs} O={bio_merged} USE_THREADING=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT", log_file)
            run_command(f"{samtools} index {bio_merged}", log_file)
            run_command(f"perl {stats_script} {bio_merged} > {dirs['merged']}/{bio_id}.bam.stats", log_file)

    print("\n✅ Pipeline complete. Log saved to:", log_file)

if __name__ == "__main__":
    main()
