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
    parser = argparse.ArgumentParser(description="Run purge_dups pipeline with ONT or PacBio reads.")
    parser.add_argument('--asmid', required=True, help='Assembly ID prefix (e.g., sample123)')
    parser.add_argument('--reads', required=True, help='ONT or PacBio reads (FASTQ or FASTQ.gz)')
    parser.add_argument('--readtype', choices=['ont', 'clr', 'hifi'], required=True, help='Type of long reads: ont, clr, or hifi')
    parser.add_argument('--threads', type=int, default=12, help='Number of threads (default: 12)')
    args = parser.parse_args()

    asm = f"{args.asmid}.fasta"
    split_fa = f"{args.asmid}.split"
    paf_file = f"{args.asmid}_reads.paf.gz"
    self_paf = f"{args.asmid}_self.paf.gz"
    cutoffs_file = f"{args.asmid}.cutoffs"
    stat_file = "PB.stat"
    cov_file = "PB.base.cov"
    dups_bed = f"{args.asmid}_dups.bed"

    # Step 1: Index the assembly
    run_cmd(f"samtools faidx {asm}")

    # Step 2: Split the assembly
    run_cmd(f"split_fa {asm} > {split_fa}")

    # Step 3: Map long reads to assembly
    if args.readtype == 'ont':
        mm2_preset = 'map-ont'
    elif args.readtype == 'clr':
        mm2_preset = 'map-pb'
    else:
        mm2_preset = 'map-hifi'

    run_cmd(
        f"minimap2 -x {mm2_preset} -t {args.threads} {asm} {args.reads} | gzip -c > {paf_file}"
    )

    # Step 4: Coverage stats and cutoff calculation
    run_cmd(f"pbcstat {paf_file}")
    run_cmd(f"calcuts {stat_file} > {cutoffs_file} 2> {args.asmid}_calcuts.log")

    # Step 5: Self-alignment of assembly
    run_cmd(f"minimap2 -xasm5 -DP -t {args.threads} {split_fa} {split_fa} | gzip -c - > {self_paf}")

    # Step 6: Run purge_dups
    run_cmd(f"purge_dups -2 -T {cutoffs_file} -c {cov_file} {self_paf} > {dups_bed} 2> {args.asmid}_purgedups.log")

    # Step 7: Extract sequences
    run_cmd(f"get_seqs -p {args.asmid} -e {dups_bed} {asm}")

    # Rename outputs from get_seqs
    os.rename("purged.fa", f"{args.asmid}.purged.fa")
    os.rename("hap.fa", f"{args.asmid}.hap.fa")

    # Step 8: Clean final output FASTAs
    for suffix in ["purged", "hap"]:
        input_fa = f"{args.asmid}.{suffix}.fa"
        output_fa = f"{args.asmid}.{suffix}.fasta"
        cmd = (
            f"awk '{{if (/^>.*/) {{print}} else {{sub(/^N*/, \"\"); sub(/N*$/, \"\"); print}}}}' {input_fa} "
            f"| seqtk seq -l0 -L100 - > {output_fa}"
        )
        run_cmd(cmd)

    print("[✔] Purge_Dups long-read pipeline complete.")


if __name__ == "__main__":
    main()
