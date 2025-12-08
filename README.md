# Evaluation of Dorado v5.2.0 *de novo* basecalling models for the detection of tRNA modifications using RNA004 chemistry

A pipeline for modification-aware base calling, alignment, and prediction of nucleotide modifications in tRNAs from ONT direct RNA sequencing reads (**SQK-RNA004**)

## Usage

### Step 1: Basecalling

Perform modification-aware sup basecalling with Dorado v5.2.0 (as of December 2025 [Dorado](https://github.com/nanoporetech/dorado) supports 8 modifications - pseU_2OmeU, m5C_2OmeC, inosine_m6A_2OmeA, 2OmeG):
```
~/dorado basecaller rna004_130bps_sup@v5.2.0 --modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG  --emit-moves --modified-bases-threshold 0 --emit-sam data/IVT.pod5 > data/IVT.sam
```
Using `--modified-bases-threshold 0` sets the reporting modification probability threshold to 0.

Convert the base called sam file to fastq with [SAMtools](https://github.com/samtools/samtools):
```
samtools fastq -T "*" data/IVT.sam > data/IVT.fastq
```

### Step 2: Alignment

For IVT tRNAs, we recommend [BWA-MEM](https://github.com/lh3/bwa) based on evaluation performed by [Lucas, M.C., Pryszcz, L.P., Medina, R. et al. 2024](https://doi.org/10.1038/s41587-023-01743-6):
```
cd ~/bwa && ./bwa mem -W 13 -k 6 -x ont2d -T 20 ~/data/reference_comb_rng.fasta ~/data/IVT.fastq > ~/data/IVT_bwa.sam
```

For *ex cellulo* tRNAs, we recommend [Parasail](https://github.com/jeffdaily/parasail), based on evaluation performed by [Sun, Y. et al. 2023](https://doi.org/10.1093/nar/gkad826):
```
~/parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f data/reference_comb_rng.fasta -q data/WT.fastq -g data/WT_parasail.sam
```

### Step 3a: Read filtering - BWA-MEM

For the IVT tRNAs aligned using BWA-MEM, filtered out reads that don't cover the entire length of tRNA to retain only full-length reads:
```
python py_scripts/filter_bwa.py data/IVT_bwa.sam
```
In this step, reads that were discarded due to mapping to more than one tRNA isoacceptor will be stored in `output/discarded_reads.csv`.

### Step 3b: Read filtering - Parasail
For the *ex cello* tRNAs that were aligned using parasail, determine the optimal alignment score threshold by:

* Aligning the reads to a reference where tRNA sequences are reversed:
```
~/parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f data/reference_comb_rng_reverse.fasta data/WT.fastq -g data/WT_parasail_reverse.sam
```

* Extracting the alignment scores of reads aligned to the reversed reference:
```
python py_scripts/extract_alignment_scores.py output/AS_thres/WT_WT-reverse_parasail_scores.csv data/WT_parasail.sam data/WT_parasail_reverse.sam
```

* Calculating the 95th percentile value of the reverse alignment score distribution:
```
python py_scripts/calc_as_thres.py output/WT_WT-reverse_parasail_scores.csv WT_parasail WT_parasail_reverse output/AS_thres/WT_thres.txt output/AS_thres/WT_thres.svg
```
This will calculate the 95 % alignment score specificity threshold.

Based on the alignment score threshold calculated in the previous step, retain the reads with alignment score above the 95 % specificity threshold, and top-scoring alignments mapping uniquely to a single tRNA isoacceptor:
```
python py_scripts/filter_parasail.py data/WT_parasail.sam 30
```
Here, 30 is the 5 % specificity alignment score threshold.
This will generate a file `data/WT_thres30_filtered.sam` with the filtered reads.
The reads that were discarded due to mapping to more than one tRNA isoacceptor will be stored in `output/discarded_reads.csv`.

Finally, retain full-length reads:
```
python py_scripts/filter_parasail_fulllen.py data/WT_thres30_filtered.sam
```
This will generate a file `data/WT_thres30_filtered_fulllength.sam` with the filtered full-length reads.

### Step 4: Dorado parsing

Assign modification score predicted by Dorado to each nucleotide in the read that aligns identically to the reference:
```
python py_scripts/parse_dorado.py data/reference_comb_rng.fasta data IVT,WT
```
This will generate a file `output/trna2mods.pckl` with the parsed results of Dorado modificiaton-calling that can be used for downstream ml-score analysis.

### Step 5: Calculating mismatch and deletion rates

Calculate the nucleotide mismatch and deletion rates at each position in reference:
```
python py_scripts/calc_mis_del_rates.py data/WT_align_final.sam data/reference_comb_rng.fasta output/mis_del_rates_WT.csv
```



