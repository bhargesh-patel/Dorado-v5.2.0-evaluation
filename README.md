# Evaluation of Dorado v5.2.0 *de novo* basecalling models for the detection of tRNA modifications using RNA004 chemistry
[![https://doi.org/10.64898/2025.12.09.693013](https://img.shields.io/badge/DOI-https://doi.org/10.64898/2025.12.09.693013-rgb%28255%2C%200%2C%200%2C%201%29)](https://doi.org/10.64898/2025.12.09.693013)

A pipeline for modification-aware basecalling, alignment, and prediction of nucleotide modifications in tRNAs from ONT direct RNA sequencing reads (**SQK-RNA004**)

## Usage

### Step 1: Basecalling

Perform modification-aware sup basecalling with Dorado v5.2.0 (as of December 2025 [Dorado](https://github.com/nanoporetech/dorado) supports eight modifications - pseU_2OmeU, m5C_2OmeC, inosine_m6A_2OmeA, 2OmeG):
```
~/dorado basecaller rna004_130bps_sup@v5.2.0 --modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG  --emit-moves --modified-bases-threshold 0 --emit-sam data/IVT.pod5 > data/IVT.sam
```
Setting `--modified-bases-threshold 0` forces Dorado to report all modification probabilities above 0.

Convert the basecalled SAM file to FASTQ using [SAMtools](https://github.com/samtools/samtools):
```
samtools fastq -T "*" data/IVT.sam > data/IVT.fastq
```

### Step 2: Alignment

For IVT tRNAs, we recommend [BWA-MEM](https://github.com/lh3/bwa), based on evaluations by [Lucas, M.C., Pryszcz, L.P., Medina, R. et al. 2024](https://doi.org/10.1038/s41587-023-01743-6):
```
cd ~/bwa && ./bwa mem -W 13 -k 6 -x ont2d -T 20 ~/data/reference_comb_rng.fasta ~/data/IVT.fastq > ~/data/IVT_bwa.sam
```

For *ex cellulo* tRNAs, we recommend [Parasail](https://github.com/jeffdaily/parasail), based on the evaluation by [Sun, Y. et al. 2023](https://doi.org/10.1093/nar/gkad826):
```
~/parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f data/reference_comb_rng.fasta -q data/WT.fastq -g data/WT_parasail.sam
```

### Step 3a: Read filtering - BWA-MEM

For IVT tRNAs aligned using BWA-MEM, filter out reads that do not cover the entire length of the tRNA, retaining only full-length reads:
```
python py_scripts/filter_bwa.py data/IVT_bwa.sam
```
Reads discarded because they map to multiple tRNA isoacceptors are stored in `output/discarded_reads.csv`.

### Step 3b: Read filtering - Parasail
For *ex cellulo* tRNAs aligned using Parasail, determine the optimal alignment-score threshold by:

* Aligning the reads to a reference where tRNA sequences are reversed:
```
~/parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f data/reference_comb_rng_reverse.fasta -q data/WT.fastq -g data/WT_parasail_reverse.sam
```

* Extracting alignment scores from reads aligned to the reversed reference:
```
python py_scripts/extract_alignment_scores.py output/AS_thres/WT_WT-reverse_parasail_scores.csv data/WT_parasail.sam data/WT_parasail_reverse.sam
```

* Calculating the 95th percentile value of the reverse alignment score distribution:
```
python py_scripts/calc_as_thres.py output/WT_WT-reverse_parasail_scores.csv WT_parasail WT_parasail_reverse output/AS_thres/WT_thres.txt output/AS_thres/WT_thres.svg
```
This calculates the 95% alignment-score specificity threshold.

Using the alignment-score threshold from the previous step, retain reads with scores above the 95% specificity threshold and uniquely mapping to a single tRNA isoacceptor:
```
python py_scripts/filter_parasail.py data/WT_parasail.sam 30
```
Here, 30 is the 95 % specificity alignment score threshold.
This generates a file `data/WT_thres30_filtered.sam` with the filtered reads.
The reads discarded due to mapping to more than one tRNA isoacceptor are stored in `output/discarded_reads.csv`.

Finally, retain full-length reads:
```
python py_scripts/filter_parasail_fulllen.py data/WT_thres30_filtered.sam
```
This generates a file `data/WT_thres30_filtered_fulllength.sam` with the filtered full-length reads.

### Step 4: Dorado parsing

After filtering and retaining full-length reads, rename the final SAM files for ease in the format `<sample>_align_final.sam`.

Assign Dorado modification scores only to nucleotides that match the reference in the final filtered full-length SAM files (`<sample>_align_final.sam`):
```
python py_scripts/parse_dorado.py data/reference_comb_rng.fasta data IVT,WT
```
This generates a file `output/trna2mods.pckl` with the parsed results of Dorado modification-calling, which can be used for downstream ml-score analysis.

### Step 5: Calculating mismatch and deletion rates

Calculate the nucleotide mismatch and deletion rates at each position in reference:
```
python py_scripts/calc_mis_del_rates.py data/WT_align_final.sam data/reference_comb_rng.fasta output/mis_del_rates_WT.csv
```



