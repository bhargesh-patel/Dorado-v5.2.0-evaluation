#!/usr/bin/env python3
import pysam
import csv
import sys
from collections import defaultdict
from Bio import SeqIO

def calculate_rates_per_trna(sam_file, fasta_file, output_csv):
    # Load reference sequences from FASTA
    reference_sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

    # Initialise counts for each tRNA
    counts = {}
    for trna_id, seq in reference_sequences.items():
        counts[trna_id] = [
            {"matches": 0, "mismatches": 0, "deletions": 0, "base": base}
            for base in seq
        ]

    # Open SAM/BAM file
    samfile = pysam.AlignmentFile(sam_file, "r")

    for read in samfile.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        trna_id = samfile.get_reference_name(read.reference_id)
        if trna_id not in counts:
            continue

        # Get aligned pairs without requiring reference sequence from SAM header
        for qpos, rpos in read.get_aligned_pairs(with_seq=False):
            if rpos is None:
                # Insertion relative to reference — skip for deletion rate
                continue
            if rpos >= len(counts[trna_id]):
                continue  # Position outside reference length

            ref_base = reference_sequences[trna_id][rpos]
            if qpos is None:
                # Deletion from query
                counts[trna_id][rpos]["deletions"] += 1
            else:
                qbase = read.query_sequence[qpos]
                if qbase.upper() == ref_base.upper():
                    counts[trna_id][rpos]["matches"] += 1
                else:
                    counts[trna_id][rpos]["mismatches"] += 1

    samfile.close()

    # Write results to CSV
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["tRNA", "position", "nucleotide", "mismatch_rate", "deletion_rate"])

        for trna_id, positions in counts.items():
            for i, pos_data in enumerate(positions):
                total = pos_data["matches"] + pos_data["mismatches"] + pos_data["deletions"]
                mismatch_rate = pos_data["mismatches"] / total if total > 0 else 0
                deletion_rate = pos_data["deletions"] / total if total > 0 else 0
                writer.writerow([
                    trna_id,
                    i + 1,  # 1-based position
                    pos_data["base"],
                    round(mismatch_rate, 4),
                    round(deletion_rate, 4)
                ])

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.sam> <reference.fasta> <output.csv>")
        sys.exit(1)

    calculate_rates_per_trna(sys.argv[1], sys.argv[2], sys.argv[3])
