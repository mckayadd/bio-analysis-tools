#!/usr/bin/env python3
"""
fasta_analyzer.py
-----------------
Analyzes DNA sequences in a multi-FASTA file:
1. Counts records
2. Finds sequence lengths, shortest/longest sequences
3. Identifies Open Reading Frames (ORFs) for a given frame
4. Finds repeats of a given length n
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
from collections import defaultdict, Counter


# ------------------------------------------------------------
# 1. Parse FASTA file and basic stats
# ------------------------------------------------------------
def parse_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    if not records:
        print("No records found. Please check your FASTA file.")
        sys.exit(1)
    return records


def fasta_summary(records):
    print(f"\n[1] Number of records in file: {len(records)}")

    lengths = {rec.id: len(rec.seq) for rec in records}
    print("\n[2] Sequence lengths:")
    for rec_id, length in lengths.items():
        print(f"  {rec_id}: {length} bp")

    max_len = max(lengths.values())
    min_len = min(lengths.values())

    longest = [rid for rid, l in lengths.items() if l == max_len]
    shortest = [rid for rid, l in lengths.items() if l == min_len]

    print(f"\n  Longest sequence(s): {', '.join(longest)} ({max_len} bp)")
    print(f"  Shortest sequence(s): {', '.join(shortest)} ({min_len} bp)")

    return lengths


# ------------------------------------------------------------
# 2. Find ORFs in a given reading frame (1, 2, or 3)
# ------------------------------------------------------------
def find_orfs(seq_str, frame=1):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    seq = seq_str[frame - 1:]  # adjust for frame
    seq_len = len(seq)

    i = 0
    while i < seq_len - 2:
        codon = seq[i:i + 3]
        if codon == start_codon:
            # found start codon
            for j in range(i + 3, seq_len - 2, 3):
                stop = seq[j:j + 3]
                if stop in stop_codons:
                    orf_seq = seq[i:j + 3]
                    start_pos = i + frame  # +frame because of offset
                    orfs.append((orf_seq, start_pos))
                    break
        i += 3
    return orfs


def longest_orf_in_records(records, frame=1):
    longest_orf = ""
    longest_id = None
    start_position = None

    for rec in records:
        orfs = find_orfs(str(rec.seq).upper(), frame)
        if orfs:
            seq_orf, start = max(orfs, key=lambda x: len(x[0]))
            if len(seq_orf) > len(longest_orf):
                longest_orf = seq_orf
                longest_id = rec.id
                start_position = start

    if longest_orf:
        print(f"\n[3] Longest ORF in reading frame {frame}:")
        print(f"  Sequence ID: {longest_id}")
        print(f"  ORF length: {len(longest_orf)} bp")
        print(f"  Starts at position: {start_position}")
        print(f"  ORF: {longest_orf}")
    else:
        print(f"\n[3] No ORFs found in reading frame {frame}.")


# ------------------------------------------------------------
# 3. Find repeats of length n
# ------------------------------------------------------------
def find_repeats(records, n):
    repeats = Counter()
    for rec in records:
        seq = str(rec.seq).upper()
        for i in range(len(seq) - n + 1):
            repeat = seq[i:i + n]
            repeats[repeat] += 1

    # Keep only those that occur more than once
    repeats = {r: c for r, c in repeats.items() if c > 1}

    if repeats:
        max_repeat = max(repeats.values())
        most_common = [r for r, c in repeats.items() if c == max_repeat]
        print(f"\n[4] Most frequent repeat(s) of length {n}:")
        for r in most_common:
            print(f"  {r} → {max_repeat} occurrences")
    else:
        print(f"\n[4] No repeats of length {n} found.")


# ------------------------------------------------------------
# Main driver
# ------------------------------------------------------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python fasta_analyzer.py <fasta_file> [frame] [repeat_length]")
        print("Example: python fasta_analyzer.py dna.example.fasta 1 6")
        sys.exit(1)

    fasta_file = sys.argv[1]
    frame = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    repeat_len = int(sys.argv[3]) if len(sys.argv) > 3 else 6

    print(f"\nAnalyzing FASTA file: {fasta_file}")

    records = parse_fasta(fasta_file)
    fasta_summary(records)
    longest_orf_in_records(records, frame)
    find_repeats(records, repeat_len)


if __name__ == "__main__":
    main()
