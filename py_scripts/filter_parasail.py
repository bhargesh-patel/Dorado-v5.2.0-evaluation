#!/usr/bin/env python3
import sys
import os
import csv
import pysam
from copy import deepcopy
from collections import defaultdict

def get_trna_identity(full_name):
    if full_name.startswith("tRNA-"):
        parts = full_name.split("-")
        if len(parts) >= 2:
            return f"{parts[0]}-{parts[1]}"
        else:
            return full_name
    elif full_name.upper().startswith("SPMITTRNA"):
        return full_name.split(".")[0].upper()
    else:
        return full_name

def filter_and_preserve_tags(input_sam, as_threshold):
    input_dir = os.path.dirname(input_sam)
    sam_basename = os.path.basename(input_sam).split("_parasail")[0]

    # Output SAM in same directory as input
    output_sam = os.path.join(input_dir, f"{sam_basename}_thres{as_threshold}_filtered.sam")

    # CSV log in output/ directory
    log_dir = "output"
    os.makedirs(log_dir, exist_ok=True)
    log_csv = os.path.join(log_dir, f"{sam_basename}_thres{as_threshold}_discarded_reads.csv")

    # Derive Dorado SAM based on input (before "_parasail")
    dorado_sam = os.path.join(input_dir, f"{sam_basename}.sam")

    samfile = pysam.AlignmentFile(input_sam, "r")
    outfile = pysam.AlignmentFile(output_sam, "w", template=samfile)
    movefile = pysam.AlignmentFile(dorado_sam, "r", check_sq=False)

    read2best = {}
    discarded_reads = {}

    # Step 1: filter based on AS
    for aln in samfile.fetch():
        try:
            as_score = aln.get_tag("AS")
        except KeyError:
            continue
        if as_score < as_threshold:
            continue

        if aln.query_name not in read2best:
            read2best[aln.query_name] = [deepcopy(aln)]
        else:
            current_max = read2best[aln.query_name][0].get_tag("AS")
            if as_score > current_max:
                read2best[aln.query_name] = [deepcopy(aln)]
            elif as_score == current_max:
                read2best[aln.query_name].append(deepcopy(aln))

    samfile.close()

    # Step 2: check for multiple tRNA identities
    for read_id, aligns in read2best.items():
        identities = set(get_trna_identity(aln.reference_name) for aln in aligns)
        if len(identities) > 1:
            discarded_reads[read_id] = set(aln.reference_name for aln in aligns)

    # Step 3: preserve tags from Dorado SAM and write kept reads
    for aln in movefile.fetch():
        if aln.query_name in read2best:
            for i in range(len(read2best[aln.query_name])):
                read2best[aln.query_name][i].set_tag("ts", aln.get_tag("ts"))
                read2best[aln.query_name][i].set_tag("mv", aln.get_tag("mv"))
                if aln.has_tag("MM"):
                    read2best[aln.query_name][i].set_tag("MM", aln.get_tag("MM"))
                    read2best[aln.query_name][i].set_tag("ML", aln.get_tag("ML"))
                if aln.has_tag("pi"):
                    read2best[aln.query_name][i].set_tag("pi", aln.get_tag("pi"))
                    read2best[aln.query_name][i].set_tag("sp", aln.get_tag("sp"))

    # Write filtered reads
    for read_id, aligns in read2best.items():
        if read_id not in discarded_reads:
            for aln in aligns:
                outfile.write(aln)

    movefile.close()
    outfile.close()

    # Write discarded reads log in output/
    with open(log_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["read_id", "aligned_tRNAs"])
        for rid, refs in discarded_reads.items():
            writer.writerow([rid, ";".join(sorted(refs))])

    print(f"Filtered SAM written to: {output_sam}")
    print(f"Discarded reads log written to: {log_csv}")
    print(f"Total reads discarded: {len(discarded_reads)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.sam> <AS_threshold>")
        sys.exit(1)

    input_sam = sys.argv[1]
    as_threshold = int(sys.argv[2])
    filter_and_preserve_tags(input_sam, as_threshold)
