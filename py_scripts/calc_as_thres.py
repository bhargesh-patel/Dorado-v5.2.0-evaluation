import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
import os

csv.field_size_limit(10**10)  # increase field size limit
def read_scores(file_path, sample1, sample2):
    delimiter = '\t' if file_path.endswith('.tsv') else ','
    scores1, scores2 = [], []

    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            sample = row['Sample']
            score_list = list(map(int, row['Scores'].split(',')))
            if sample == sample1:
                scores1 = score_list
            elif sample == sample2:
                scores2 = score_list

    if not scores1 or not scores2:
        print("Error: One or both sample names not found in the file.")
        sys.exit(1)

    return scores1, scores2

def compute_95_spec_threshold(real_scores, random_scores):
    real_scores = np.array(real_scores)
    random_scores = np.array(random_scores)

    sorted_random = np.sort(random_scores)
    index = int(0.95 * len(sorted_random)) - 1
    index = max(0, min(index, len(sorted_random) - 1))
    threshold = sorted_random[index]
    return threshold

def plot_histograms(real_scores, random_scores, threshold, output_path):
    plt.figure(figsize=(10, 6))
    bins = np.linspace(min(real_scores + random_scores), max(real_scores + random_scores), 100)

    plt.hist(real_scores, bins=bins, alpha=0.5, label='Real (Forward)', color='skyblue', density=True)
    plt.hist(random_scores, bins=bins, alpha=0.5, label='Random (Reverse)', color='salmon', density=True)

    plt.axvline(threshold, color='black', linestyle='--', linewidth=2, label=f'Threshold = {threshold}')
    plt.xlabel("Alignment Score")
    plt.ylabel("Density")
    plt.title("Alignment Score Distribution with 95% Specificity Threshold")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, format='svg')
    plt.close()

def main():
    if len(sys.argv) != 6:
        print("Usage: python calc_threshold_plot.py scores.csv sample1 sample2 output_threshold.txt output_plot.svg")
        sys.exit(1)

    score_file = sys.argv[1]
    sample1 = sys.argv[2]
    sample2 = sys.argv[3]
    threshold_output_file = sys.argv[4]
    plot_output_file = sys.argv[5]

    real_scores, random_scores = read_scores(score_file, sample1, sample2)
    threshold = compute_95_spec_threshold(real_scores, random_scores)

    with open(threshold_output_file, 'w') as f:
        f.write(f"{threshold}\n")

    print(f"✅ Threshold at 95% specificity: {threshold}")
    print(f"📄 Written to: {threshold_output_file}")

    plot_histograms(real_scores, random_scores, threshold, plot_output_file)
    print(f"📈 Histogram saved to: {plot_output_file}")

if __name__ == "__main__":
    main()
