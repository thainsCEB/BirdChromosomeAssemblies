#!/usr/bin/env python3

import os
import argparse
import subprocess
import pandas as pd

def run_cmd(cmd, log_file):
    with open(log_file, 'a') as log:
        log.write(f"Running command: {' '.join(cmd)}\n")
        subprocess.run(cmd, stdout=log, stderr=log, check=True)

def align_pair(read1, read2, reference, threads, outbam, log):
    bwa_cmd = [
        "bwa", "mem", "-SP5M", "-t", str(threads),
        reference, read1, read2
    ]
    view_cmd = ["samtools", "view", "-bSh", "-"]
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", outbam]

    bwa = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
    view = subprocess.Popen(view_cmd, stdin=bwa.stdout, stdout=subprocess.PIPE)
    sort = subprocess.Popen(sort_cmd, stdin=view.stdout)
    sort.communicate()

def map_hic_reads(args):
    os.makedirs(args.outdir, exist_ok=True)
    log_path = os.path.join(args.outdir, "mapping.log")
    df = pd.read_csv(args.map, sep="\t", header=None,
                     names=["read1", "read2", "replicate_id", "sample_id"])

    replicate_groups = df.groupby("replicate_id")
    for rep_id, group in replicate_groups:
        sample_ids = group["sample_id"].unique()
        if len(sample_ids) > 1:
            raise ValueError(f"Replicate ID {rep_id} maps to multiple sample IDs: {sample_ids}")
        sample_id = sample_ids[0]

        bam_files = []
        print(f"\nProcessing replicate group: {rep_id} for sample: {sample_id}")
        for idx, row in group.iterrows():
            read1, read2 = row.read1, row.read2
            align_tag = f"{sample_id}_{rep_id}_part{idx}"
            bam_path = os.path.join(args.outdir, f"{align_tag}.bam")
            print(f"Aligning: {read1} & {read2}")
            align_pair(read1, read2, args.reference, args.threads, bam_path, log_path)
            bam_files.append(bam_path)

        # Merge BAMs
        merged_bam = os.path.join(args.outdir, f"{sample_id}_{rep_id}.merged.bam")
        merge_cmd = ["samtools", "merge", "-@", str(args.threads), merged_bam] + bam_files
        run_cmd(merge_cmd, log_path)

        # Deduplication
        dedup_bam = os.path.join(args.outdir, f"{sample_id}_{rep_id}.dedup.bam")
        metrics_file = os.path.join(args.outdir, f"{sample_id}_{rep_id}.metrics.txt")
        picard_cmd = [
            "picard", "MarkDuplicates",
            f"I={merged_bam}",
            f"O={dedup_bam}",
            f"M={metrics_file}",
            "REMOVE_DUPLICATES=true",
            "VALIDATION_STRINGENCY=LENIENT"
        ]
        run_cmd(picard_cmd, log_path)

        # Index
        index_cmd = ["samtools", "index", dedup_bam]
        run_cmd(index_cmd, log_path)

        # Flagstat
        stats_file = os.path.join(args.outdir, f"{sample_id}_{rep_id}.flagstat.txt")
        flagstat_cmd = ["samtools", "flagstat", dedup_bam]
        with open(stats_file, 'w') as stats:
            subprocess.run(flagstat_cmd, stdout=stats, check=True)

        print(f"Finished sample {sample_id} (replicate {rep_id})\n")

def main():
    parser = argparse.ArgumentParser(description="Hi-C Mapping with Replicate Merging and Sample Labels")
    parser.add_argument("--map", required=True, help="Tab-separated file: read1, read2, replicate_id, sample_id")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA (BWA-indexed)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()
    map_hic_reads(args)

if __name__ == "__main__":
    main()
