import sys
import csv
import os

def extract_alignment_scores(sam_file, included_refs):
    scores = []
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            rname = fields[2]
            if any(skip in rname for skip in ["intron", "Gln", "SPMIT"]):
                continue
            for field in fields:
                if field.startswith('AS:i:'):
                    try:
                        score = int(field.split(':')[2])
                        scores.append(score)
                        included_refs.add(rname)
                    except ValueError:
                        continue
    return scores

def main():
    if len(sys.argv) < 3:
        print("Usage: python extract_alignment_scores.py output.csv input1.sam input2.sam ...")
        sys.exit(1)

    output_file = sys.argv[1]
    sam_files = sys.argv[2:]

    if output_file.endswith(".tsv"):
        delimiter = '\t'
    elif output_file.endswith(".csv"):
        delimiter = ','
    else:
        print("Error: Output file must end with .csv or .tsv")
        sys.exit(1)

    included_refs = set()
    rows = []

    for sam in sam_files:
        sample_name = os.path.splitext(os.path.basename(sam))[0]
        scores = extract_alignment_scores(sam, included_refs)
        score_str = ",".join(map(str, scores))
        rows.append([sample_name, score_str])

    # Write R-friendly output
    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter=delimiter)
        writer.writerow(["Sample", "Scores"])
        writer.writerows(rows)

    # Print unique references used
    print("✅ Written:", output_file)
    print("📋 Unique tRNAs with scores:")
    for r in sorted(included_refs):
        print(" -", r)

if __name__ == "__main__":
    main()
