#!/usr/bin/env python3
import pysam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pod5
import argparse
import pickle
import os
from scipy.stats import ks_2samp

# Adapter sequences (5' adapter is 23 nt, as in original script)
ADAPT_5P = "CCTAAGAGCAAGAAGAAGCCTGG"
ADAPT_3P = "GGCTTCTTCTTGCTCTTAGG"
ADAPT_5P_LEN = len(ADAPT_5P)  # 23

def get_signal(dim, x, signal_calibrated, signal_cpts, stride):
    """Extract per-base signal aligned to reference from SAM move table."""
    signal = np.zeros(dim)
    seq2ref = dict(x.get_aligned_pairs())

    for base_i in range(1, x.query_length):
        if seq2ref[base_i] is not None:
            ref_pos = seq2ref[base_i]
            chunk = signal_calibrated[signal_cpts[-base_i-1]:signal_cpts[-base_i]]
            step_size = max(1, len(chunk) // stride)
            for step_i in range(stride):
                subchunk = chunk[step_i*step_size:(step_i+1)*step_size]
                signal[ref_pos * stride + step_i] = np.median(subchunk)
    return signal

def extract_signal(args):
    """Extract signals per tRNA from SAM and POD5 files and save as pickle."""
    reference = pysam.FastaFile(args.ref)
    samples = args.samples.split(",")
    pod5s = args.pod5_dir.split(",")
    n_sub = int(args.subsample)

    ref2sigs = {}
    for i, sample in enumerate(samples):
        pod5_file = pod5s[i]
        print(f"Processing sample: {sample}")
        sam_path = os.path.join(args.sam_dir, sample + "_align_final.sam")
        samfile = pysam.AlignmentFile(sam_path, "r")
        pod5file = pod5.DatasetReader(pod5_file)

        ref2sigs[sample] = {ref: [] for ref in reference.references}

        for x in samfile.fetch():
            if args.tRNA and x.reference_name != args.tRNA:
                continue
            if len(ref2sigs[sample][x.reference_name]) >= n_sub:
                continue

            stride, move, ts = x.get_tag("mv")[0], np.array(x.get_tag("mv")[1:]), x.get_tag("ts")
            move_cpts = np.where(move == 1)[0]
            signal_cpts = ts + move_cpts * stride

            read_record = pod5file.get_read(x.query_name)
            if read_record is not None:
                signal_calibrated = read_record.calibrate_signal_array(read_record.signal)
            else:
                parent_id = x.get_tag("pi")
                parent_record = pod5file.get_read(parent_id)
                signal_calibrated = parent_record.calibrate_signal_array(parent_record.signal)
                signal_cpts += x.get_tag("sp")

            ref_seq = reference.fetch(x.reference_name)
            ref_signal = get_signal(stride * len(ref_seq), x, signal_calibrated, signal_cpts, stride)
            matches = set(list(zip(*x.get_aligned_pairs(matches_only=True)))[1])
            ref2sigs[sample][x.reference_name].append((ref_signal, matches))

        samfile.close()

    os.makedirs("output", exist_ok=True)
    out_file = f"output/signals_{args.samples}_{args.subsample}.pckl"
    with open(out_file, 'wb') as handle:
        pickle.dump(ref2sigs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved pickle file to {out_file}")

def get_signal_median(sample, trna_ref, ref2sigs):
    """Compute median signal per position for a sample/trna."""
    signals = [x[0] for x in ref2sigs[sample][trna_ref] if np.any(x[0] != 0)]
    if not signals:
        raise ValueError(f"No valid signals found for {trna_ref} in sample {sample}")
    signal_array = np.vstack(signals)
    signal_array_masked = np.ma.masked_where(signal_array == 0, signal_array)
    return np.ma.median(signal_array_masked, axis=0).filled(0)

def plot_signal(args, stride=6):
    """Plot signal matching original script's Excel handling, indexing, y-bounds, and legend."""
    wt, mt = args.sample1, args.sample2
    print("WT:", wt, "MT:", mt)

    with open(args.signals, 'rb') as handle:
        ref2sigs = pickle.load(handle)

    # Map Excel tRNA name
    trna_excel = args.trna
    if trna_excel.startswith("tRNA-"):
        trna_excel = trna_excel[5:]
    if trna_excel.endswith("-1"):
        trna_excel = trna_excel[:-2]
    trna_ref = args.trna

    kmer_len = args.kmer
    half_k = kmer_len // 2
    pos_of_interest = args.pos

    # Read Excel and filter by tRNA
    df_annotation = pd.read_excel(args.annotation)
    df_trna = df_annotation[df_annotation["tRNA"] == trna_excel]
    
    sel_row = df_trna[df_trna["position"] == pos_of_interest]
    if sel_row.empty:
        raise ValueError(f"Position {pos_of_interest} not found for tRNA {trna_excel} in Excel annotation")
    
    # Drop NA rows (matching original)
    df_trna_cleared = df_trna.dropna(subset=["nucleotide"], ignore_index=True)
    row_pos = df_trna_cleared[df_trna_cleared["position"] == pos_of_interest].index[0]

    # Construct k-mer sequence with adapters
    kmer_seq = ""
    for i in range(-half_k, half_k + 1):
        ref_idx = row_pos + i
        if ref_idx < 0:
            kmer_seq += ADAPT_5P[ADAPT_5P_LEN + ref_idx]
        elif ref_idx >= len(df_trna_cleared):
            kmer_seq += ADAPT_3P[ref_idx - len(df_trna_cleared)]
        else:
            kmer_seq += df_trna_cleared.iloc[ref_idx]["nucleotide"]
    print("k-mer sequence:", kmer_seq)

    # Compute offset
    diff_median = get_signal_median(mt, trna_ref, ref2sigs) - get_signal_median(wt, trna_ref, ref2sigs)
    offset = np.median(diff_median)
    print("offset =", offset)
    
    # Plot styling (matching original)
    sns.set(rc={'figure.figsize': (10, 6)})   
    sns.set_context("paper", font_scale=2)
    sns.set_style("ticks")

    # Reference indexing (matching original: 23 + row_pos)
    kmer_i = ADAPT_5P_LEN + row_pos
    length = kmer_len * stride
    x, y, hue = [], [], []
    
    print(f"Excel position: {pos_of_interest}")
    print(f"DataFrame row_pos: {row_pos}")
    print(f"Signal kmer_i: {kmer_i}")

    for sample in [wt, mt]:
        for sig, matches in ref2sigs[sample][trna_ref]:
            # Exact original slice indexing
            start = (kmer_i - half_k - 2) * stride
            end = (kmer_i + half_k - 1) * stride
            signal_chunk = sig[start:end]

            if 0 not in signal_chunk:
                if sample == wt:
                    y.extend(list(signal_chunk + offset))
                else:
                    y.extend(list(signal_chunk))
                x.extend(list(np.arange(length)))
                hue.extend([sample] * length)

    if not y:
        raise ValueError("No valid signal chunks to plot.")

    df = pd.DataFrame({"x": x, "y": y, "sample": hue})
    g = df.groupby(["sample", "x"])["y"]
    summary = g.agg(["mean", "std"]).reset_index()
    summary["std"] = summary["std"].fillna(0.0)
    summary["lower"] = summary["mean"] - summary["std"]
    summary["upper"] = summary["mean"] + summary["std"]
    shade_min = float(summary["lower"].min())
    shade_max = float(summary["upper"].max())

    y_bottom = 0
    if shade_max > y_bottom:
        pad = 0.1 * (shade_max - y_bottom)  # 10% headroom
    else:
        pad = 10
    y_top = shade_max + pad

    # Plot (matching original)
    sns.lineplot(x=x, y=y, hue=hue, errorbar="sd",
                 palette=["#008080","#FF7F50"], lw=2, alpha=0.7)
    plt.vlines(np.arange(stride, length, stride), ymin=0, ymax=y_top,
               ls="--", color="grey", lw=0.5)
    
    for base_i in range(len(kmer_seq)):
        plt.text(base_i * stride + 1, y_bottom + 0.85 * (y_top - y_bottom),
                 kmer_seq[base_i])

    plt.ylabel("Signal intensity")
    plt.xticks([])
    
    # y-limits and ticks based on shaded area
    plt.ylim(y_bottom, y_top)
    y_ticks = np.arange(y_bottom, y_top + 1e-9, 25)
    plt.yticks(y_ticks)

    sns.despine()

    # Save (matching original filename logic)
    os.makedirs("output/Signal_plots", exist_ok=True)
    if trna_ref.startswith("SPMIT") or trna_ref.endswith("intron"):
        fname = f"{wt}_{mt}_{trna_ref}_{pos_of_interest}_{kmer_len}mer.svg"
    else:
        fname = f"{wt}_{mt}_{trna_ref[5:-2]}_{pos_of_interest}_{kmer_len}mer.svg"
    out_path = os.path.join("output/Signal_plots", fname)
    plt.savefig(out_path, format="svg", bbox_inches='tight')
    print(f"Saved plot to {out_path}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Extract
    extract_parser = subparsers.add_parser("extract", help="Extract signals from a list of samples")
    extract_parser.add_argument("--samples", type=str, required=True, help="Comma-separated list of sample names")
    extract_parser.add_argument("--ref", type=str, required=True, help="Reference FASTA file")
    extract_parser.add_argument("--pod5_dir", type=str, required=True, help="Comma-separated list of POD5 file paths")
    extract_parser.add_argument("--subsample", type=str, required=True, help="Maximum number of reads per tRNA")
    extract_parser.add_argument("--sam_dir", type=str, required=True, help="Directory containing SAM files")
    extract_parser.add_argument("--tRNA", type=str, default=None, help="Specific tRNA to extract (optional)")
    extract_parser.add_argument("--annotation", type=str, required=True, help="Excel annotation file")
    extract_parser.set_defaults(func=extract_signal)

    # Plot
    plot_parser = subparsers.add_parser("plot", help="Plot signals for a specific tRNA position")
    plot_parser.add_argument("--sample1", type=str, required=True, help="Wildtype sample")
    plot_parser.add_argument("--sample2", type=str, required=True, help="Mutant sample")
    plot_parser.add_argument("--signals", type=str, required=True, help="Pickle file with signals")
    plot_parser.add_argument("--trna", type=str, required=True, help="tRNA name (e.g., tRNA-Gln-TTG-3-1)")
    plot_parser.add_argument("--pos", type=int, required=True, help="Position from Excel")
    plot_parser.add_argument("--kmer", type=int, required=True, help="K-mer length")
    plot_parser.add_argument("--annotation", type=str, required=True, help="Excel annotation file")
    plot_parser.set_defaults(func=plot_signal)

    args = parser.parse_args()
    args.func(args)
