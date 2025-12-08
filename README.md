# Evaluation of Dorado v5.2.0 *de novo* basecalling models for the detection of tRNA modifications using RNA004 chemistry

A pipeline for modification-aware base calling, alignment, and prediction of nucleotide modifications in tRNAs from ONT direct RNA sequencing reads (**SQK-RNA004**)

## Usage

### Step 1: Basecalling

To perform modification-aware sup basecalling with Dorado v5.2.0 (as of December 2025 [Dorado](https://github.com/nanoporetech/dorado) supports 8 modifications - pseU_2OmeU, m5C_2OmeC, inosine_m6A_2OmeA, 2OmeG):
```
~/dorado basecaller rna004_130bps_sup@v5.2.0 --modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG  --emit-moves --modified-bases-threshold 0 --emit-sam data/IVT.pod5 > data/IVT.sam
```
`--modified-bases-threshold 0` sets the reporting modification probability threshold to 0.

### Step 2: Convert the base called sam file to fastq with [SAMtools](https://github.com/samtools/samtools):
```
samtools fastq -T "*" —- data/IVT.sam > data/IVT.fastq
'''

### Step 3: Alignment

For IVT tRNAs, we recommend [BWA-MEM](https://github.com/lh3/bwa) based on evaluation performed by [Lucas, M.C., Pryszcz, L.P., Medina, R. et al. 2024](https://doi.org/10.1038/s41587-023-01743-6):
```
cd ~/bwa && ./bwa mem -W 13 -k 6 -x ont2d -T 20 ~/data/reference_comb_rng.fasta ~/data/IVT.fastq > ~/data/IVT_bwa.sam
```

For *ex cellulo* tRNAs, we recommend [Parasail](https://github.com/jeffdaily/parasail), based on evaluation performed by [Sun, Y. et al. 2023](https://doi.org/10.1093/nar/gkad826):
```
~/parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f data/reference_comb_rng.fasta data/ 
```
