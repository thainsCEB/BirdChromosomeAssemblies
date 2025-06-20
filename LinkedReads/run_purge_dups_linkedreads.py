#!/usr/bin/env python3

import os
import subprocess
import argparse


def run_cmd(cmd, logfile=None):
    print(f"[Running] {cmd}")
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True) as process:
        with open(logfile, 'w') if logfile else open(os.devnull, 'w') as logf:
            for line in process.stdout:
                print(line, end='')
                logf.write(line)
        process.wait()
        if process.returncode != 0:
            raise RuntimeError(f"Command failed: {cmd}")


def main():
    parser = argparse.ArgumentParser(description="Run purge_dups pipeline with linked-reads.")
    parser.add_argument('--asmid', required=True, help='Assembly ID prefix (e.g., sample123)')
    parser.add_argument('--r1', required=True, help='Paired R1 fastq.gz')
    parser.add_argument('--r2', required=True, help='Paired R2 fastq.gz')
    parser.add_argument('--unpaired', help='Unpaired reads (optional)')
    parser.add_argument('--threads', type=int, default=12, help='Number of threads (default: 12)')
    args = parser.parse_args()

    asm = f"{args.asmid}_supernova.fasta"

    # Step 1: Index the reference
    run_cmd(f"bwa index {asm}")
    run_cmd(f"samtools faidx {asm}")

    # Step 2: Align reads
    paired_bam = f"{args.asmid}_purgeDups_paired.bam"
    run_cmd(f"bwa mem -t {args.threads} {asm} {args.r1} {args.r2} | samtools view -@ {args.threads} -bSh - > {paired_bam}")

    if args.unpaired:
        single_bam = f"{args.asmid}_purgeDups_single.bam"
        run_cmd(f"bwa mem -t {args.threads} {asm} {args.unpaired} | samtools view -@ {args.threads} -hbS - > {single_bam}")
        merged_bam = f"{args.asmid}_purgeDups.bam"
        run_cmd(f"samtools merge -@ {args.threads} {merged_bam} {paired_bam} {single_bam}")
    else:
        merged_bam = paired_bam

    # Step 3: Calculate statistics
    run_cmd(f"ngscstat {merged_bam}")
    run_cmd(f"calcuts TX.stat > {args.asmid}.cutoffs 2> {args.asmid}_calcuts.log")

    # Step 4: Split FASTA
    split_fa = f"{args.asmid}_supernova.split"
    run_cmd(f"split_fa {asm} > {split_fa}")

    # Step 5: Align assembly to itself
    paf = f"{args.asmid}_supernova.split.paf.gz"
    run_cmd(f"minimap2 -xasm5 -DP {split_fa} {split_fa} | gzip -c - > {paf}")

    # Step 6: Purge dups
    run_cmd(f"purge_dups -2 -T {args.asmid}.cutoffs -c TX.base.cov {paf} > {args.asmid}_dups.bed 2> {args.asmid}_purgedups.log")

    # Step 7: Extract purged sequences
    run_cmd(f"get_seqs -p {args.asmid} -e {args.asmid}_dups.bed {asm}")

    # Step 8: Clean and filter output FASTAs
    for suffix in ["purged", "hap"]:
        input_fa = f"{args.asmid}.{suffix}.fa"
        output_fa = f"{args.asmid}.supernova_{suffix}.fasta"
        cmd = (
            f"awk '{{if (/^>.*/) {{print}} else {{sub(/^N*/, \"\"); sub(/N*$/, \"\"); print}}}}' {input_fa} "
            f"| seqtk seq -l0 -L100 - > {output_fa}"
        )
        run_cmd(cmd)

    print("[✔] Purge_Dups pipeline complete.")


if __name__ == "__main__":
    main()
