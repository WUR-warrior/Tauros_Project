#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import argparse
import random
import numpy as np
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare haplotypes with background similarity correction.\n"
        )
    )
    parser.add_argument('--vcf', required=True, help='Phased VCF file (.vcf.gz)')
    parser.add_argument('--admixed', required=True, help='Sample ID of the admixed individual')
    parser.add_argument('--references', required=True, nargs='+', help='List of reference sample IDs')
    parser.add_argument('--window_size', type=int, default=100,
                        help='Window size in SNPs (suggested: 100-500 based on LD structure)')
    parser.add_argument('--background_samples', type=int, default=1000, help='Number of background samples for null distribution')
    parser.add_argument('--z_threshold', type=float, default=1.0, help='Z-score threshold for reporting match')
    parser.add_argument('--out', required=True, help='Path to output CSV file')
    return parser.parse_args()

def get_sample_indices(header_line, target_samples):
    samples = header_line.strip().split('\t')[9:]
    indices = {s: i for i, s in enumerate(samples) if s in target_samples}
    return indices, samples

def extract_genotypes(vcf_path, target_indices):
    haplotypes = []
    positions = []

    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                break
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos = fields[0], fields[1]
            genotypes = [fields[9 + i].split(':')[0] for i in range(len(target_indices))]
            alleles = [g.replace('|', '/') for g in genotypes]
            if all('/' in g for g in alleles):
                haplotypes.append([a.split('/') for a in alleles])
                positions.append((chrom, pos))
    return haplotypes, positions

def get_background_distribution(haplotypes, ref_index, admix_index, window_size, n_samples):
    match_counts = []
    total_snps = len(haplotypes)
    for _ in range(n_samples):
        start = random.randint(0, total_snps - window_size)
        matches = 0
        for i in range(start, start + window_size):
            a0, a1 = haplotypes[i][admix_index]
            r0, r1 = haplotypes[i][ref_index]
            if a0 == r0 or a0 == r1:
                matches += 1
            if a1 == r0 or a1 == r1:
                matches += 1
        match_counts.append(matches / (2 * window_size))
    return np.mean(match_counts), np.std(match_counts)

def main():
    args = parse_args()
    all_samples = set(args.references + [args.admixed])

    with gzip.open(args.vcf, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                sample_indices, all_sample_names = get_sample_indices(line, all_samples)
                break

    if args.admixed not in sample_indices:
        raise ValueError(f"Admixed sample '{args.admixed}' not found in VCF header.")

    admix_idx = list(sample_indices).index(args.admixed)
    ref_names = [s for s in sample_indices if s != args.admixed]

    print("Extracting phased genotypes...")
    haplotypes, positions = extract_genotypes(args.vcf, sample_indices)

    print("Computing background match distributions per reference sample...")
    background = {}
    for ref in ref_names:
        ref_idx = list(sample_indices).index(ref)
        mu, sigma = get_background_distribution(haplotypes, ref_idx, admix_idx, args.window_size, args.background_samples)
        background[ref] = (mu, sigma)

    print("Sliding through haplotype windows...")
    results = []
    num_windows = len(haplotypes) // args.window_size

    for w in range(num_windows):
        start = w * args.window_size
        end = start + args.window_size
        region = f"{positions[start][0]}:{positions[start][1]}-{positions[end-1][1]}"
        window_data = {'Region': region}

        for ref in ref_names:
            ref_idx = list(sample_indices).index(ref)
            matches = 0
            for i in range(start, end):
                a0, a1 = haplotypes[i][admix_idx]
                r0, r1 = haplotypes[i][ref_idx]
                if a0 == r0 or a0 == r1:
                    matches += 1
                if a1 == r0 or a1 == r1:
                    matches += 1
            match_fraction = matches / (2 * args.window_size)
            mu, sigma = background[ref]
            z = (match_fraction - mu) / sigma if sigma > 0 else 0
            window_data[ref] = round(z, 2) if z >= args.z_threshold else 'NA'

        results.append(window_data)

    print(f"Writing Z-score haplotype matches to: {args.out}")
    with open(args.out, 'w', newline='') as out_csv:
        writer = csv.DictWriter(out_csv, fieldnames=['Region'] + ref_names, delimiter=',')
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print("Done.")

if __name__ == '__main__':
    main()
