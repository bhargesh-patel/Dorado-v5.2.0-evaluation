#!/usr/bin/env python3
import pysam
import sys
import os
import csv
from collections import defaultdict

def get_trna_identity(trna_name):
    """
    Extract tRNA identity from full name:
    - tRNA-Ala-AGC-1-1 → tRNA-Ala
    - tRNA-Ala-AGC-1-1_intron → tRNA-Ala
    - SPMITTRNAHIS.01.1 → SPMITTRNAHIS
    """
    if trna_name.startswith("tRNA-"):
        return "-".join(trna_name.split("-")[:2])
    elif trna_name.startswith("SPMITTRNA"):
        return trna_name.split(".")[0]
    return trna_name  # fallback

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <input.sam>")
    sys.exit(1)

input_file = sys.argv[1]
input_dir = os.path.dirname(input_file)
basename = os.path.splitext(os.path.basename(input_file))[0]
output_sam = os.path.join(input_dir, f"{basename}_fulllen.sam")
output_csv = os.path.join("output", f"{basename}_discarded_reads.csv")

# Make sure output dir exists for CSV
os.makedirs("output", exist_ok=True)

dorado_file = input_file.split("_bwa_comb_rng.sam")[0] + ".sam"

samfile = pysam.AlignmentFile(input_file, "r")
outfile = pysam.AlignmentFile(output_sam, 'w', template=samfile)
movefile = pysam.AlignmentFile(dorado_file, "r", check_sq=False)

# Map Dorado reads for tag transfer
read2move = {x.query_name: x for x in movefile.fetch()}

reads_dict = defaultdict(list)
for read in samfile.fetch():
    if not read.is_unmapped:
        reads_dict[read.query_name].append(read)

discarded_reads = {}

for read_name, alignments in reads_dict.items():
    # Filter alignments based on full-length criteria
    filtered_aligns = []
    full_length_alignments = []
    for read in alignments:
        ref_pos = read.reference_start
        cigar = read.cigartuples
        first_match_pos = None
        last_match_pos = ref_pos

        for op, length in cigar:
            if op in [0, 7, 8]:  # M, =, X
                if first_match_pos is None:
                    first_match_pos = ref_pos
                ref_pos += length
                last_match_pos = ref_pos
            elif op in [2, 3]:  # D, N
                ref_pos += length
            elif op in [1, 4, 5]:  # I, S, H
                continue
            else:
                continue

        reference_length = samfile.get_reference_length(read.reference_name)

        if first_match_pos is None:
            continue

        # Keep if full length-ish
        if first_match_pos < 25 and (reference_length - last_match_pos) <= 20:
            full_length_alignments.append(read)

    if not full_length_alignments:
        continue  # no full length alignments, skip

    if len(full_length_alignments) == 1:
        filtered_aligns = full_length_alignments
    else:
        # Multiple full-length alignments: check tRNA identity uniqueness
        identities = {get_trna_identity(aln.reference_name) for aln in full_length_alignments}
        full_names = {aln.reference_name for aln in full_length_alignments}
        if len(identities) == 1:
            filtered_aligns = full_length_alignments
        else:
            discarded_reads[read_name] = full_names
            continue  # discard all alignments for this read

    # Write filtered alignments and transfer tags
    for read in filtered_aligns:
        if read.query_name in read2move:
            move_read = read2move[read.query_name]
            for tag in ["ts", "mv", "MM", "ML", "pi", "sp", "ns"]:
                if move_read.has_tag(tag):
                    read.set_tag(tag, move_read.get_tag(tag))
        outfile.write(read)

samfile.close()
movefile.close()
outfile.close()

# Write discarded reads CSV log
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["read_id", "aligned_tRNAs"])
    for rid, tnames in discarded_reads.items():
        writer.writerow([rid, ";".join(sorted(tnames))])

print(f"Filtered SAM written to: {output_sam}")
print(f"Discarded reads CSV written to: {output_csv}")
print(f"Total discarded reads: {len(discarded_reads)}")
