import pickle
import pysam
import sys

reference = pysam.FastaFile(sys.argv[1])
align_dir = sys.argv[2]
samples = sys.argv[3].split(",")

def filter_mm_ml_for_hardclipping(aln):
    if not aln.has_tag("MM") or not aln.has_tag("ML"):
        return {}, {}

    mm_tag = aln.get_tag("MM")
    ml_tag = aln.get_tag("ML")

    cigar = aln.cigartuples  # list of (operation, length)
    seq = aln.seq
    retained_base_indices = {"A": [], "C": [], "G": [], "T": []}
    base_index = {"A": 0, "C": 0, "G": 0, "T": 0}
    read_pos = 0

    for op, length in cigar:
        if op in (0, 4, 7, 8):  # M, S, =, X
            for i in range(length):
                if read_pos >= len(seq):
                    break
                base = seq[read_pos]
                if base in retained_base_indices:
                    retained_base_indices[base].append(base_index[base])
                    base_index[base] += 1
                read_pos += 1
        elif op == 1:  # Insertion
            read_pos += length
        elif op in (2, 3, 5):  # D, N, H
            continue

    index_dict, ML_dict = {}, {}
    chunk = 0
    for mod_entry in mm_tag.strip(";").split(";"):
        base = mod_entry[0]
        if base not in retained_base_indices:
            continue
        mod_key = mod_entry.split(",")[0]
        skips = list(map(int, mod_entry.split(",")[1:]))
        loc = -1
        retained_indices = []
        retained_ml = []

        for skip in skips:
            loc += skip + 1
            if loc in retained_base_indices[base]:
                retained_indices.append(loc)
                if chunk + len(retained_ml) < len(ml_tag):
                    retained_ml.append(ml_tag[chunk + len(retained_ml)])

        if retained_indices:
            new_skips = []
            last = -1
            for idx in retained_indices:
                new_skips.append(idx - last - 1)
                last = idx
            index_dict[mod_key] = retained_indices
            ML_dict[mod_key] = retained_ml

        chunk += len(skips)

    return index_dict, ML_dict

trna2mods, trna2cov = {}, {}
for sample in samples:
    print(sample)
    samfile = pysam.AlignmentFile(f"{align_dir}/{sample}_align_final.sam", "r")
    trna2mods[sample], trna2cov[sample] = {}, {}
    for ref in reference.references:
        trna2mods[sample][ref] = {}
        trna2cov[sample][ref] = 0
        ref_seq = reference.fetch(ref)
        base2pos = {
            "T+17802.": [pos for pos, char in enumerate(ref_seq) if char == "T"],
            "A+17596.": [pos for pos, char in enumerate(ref_seq) if char == "A"],
            "A+a.": [pos for pos, char in enumerate(ref_seq) if char == "A"],
            "C+m.": [pos for pos, char in enumerate(ref_seq) if char == "C"],
            "A+69426.": [pos for pos, char in enumerate(ref_seq) if char == "A"],
            "C+19228.": [pos for pos, char in enumerate(ref_seq) if char == "C"],
            "T+19227.": [pos for pos, char in enumerate(ref_seq) if char == "T"],
            "G+19229.": [pos for pos, char in enumerate(ref_seq) if char == "G"]
        }
        for mod in base2pos:
            trna2mods[sample][ref][mod] = {pos: [] for pos in base2pos[mod]}

    count = 0
    for x in samfile.fetch():
        trna2cov[sample][x.reference_name] += 1
        count += 1
        if count % 500000 == 0:
            print(count)

        try:
            index_dict, ML_dict = filter_mm_ml_for_hardclipping(x)
        except Exception as e:
            print(f"Failed MM/ML parsing for {x.query_name}: {e}")
            continue

        seq2ref = dict(x.get_aligned_pairs(matches_only=False))
        base2pos_seq = {
            "T": [pos for pos, char in enumerate(x.seq) if char == "T"],
            "A": [pos for pos, char in enumerate(x.seq) if char == "A"],
            "C": [pos for pos, char in enumerate(x.seq) if char == "C"],
            "G": [pos for pos, char in enumerate(x.seq) if char == "G"]
        }

        for mod in base2pos:
            base = mod[0]
            mod_scores = {}
            if mod in index_dict:
                for i, idx in enumerate(index_dict[mod]):
                    if idx >= len(base2pos_seq[base]):
                        continue
                    seq_pos = base2pos_seq[base][idx]
                    ref_pos = seq2ref.get(seq_pos)
                    if ref_pos is not None:
                        mod_scores[ref_pos] = ML_dict[mod][i]

            for i, seq_pos in enumerate(base2pos_seq[base]):
                ref_pos = seq2ref.get(seq_pos)
                if ref_pos is not None and reference.fetch(x.reference_name, ref_pos, ref_pos + 1) == base:
                    score = mod_scores.get(ref_pos, 0)
                    trna2mods[sample][x.reference_name][mod][ref_pos].append(score)

    samfile.close()

with open("./output/trna2mods.pckl", 'wb') as pfile:
    pickle.dump(trna2mods, pfile, protocol=pickle.HIGHEST_PROTOCOL)
