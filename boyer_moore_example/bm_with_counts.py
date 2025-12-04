#!/usr/bin/env python

"""
bm_with_counts.py

Naive exact matching and Boyer–Moore string matching
instruments that also count:
  (a) number of character comparisons
  (b) number of alignments tried
"""

from bm_preproc import BoyerMoore
from kmer_index import Index  # assumes the given code is in kmer_index.py

def naive_with_counts(p, t):
    """
    Naive exact matching.

    Args:
        p: pattern (string)
        t: text (string)

    Returns:
        occurrences: list of starting indices where p matches t
        num_alignments: how many alignments of p against t were tried
        num_comparisons: how many character comparisons were performed
    """
    occurrences = []
    num_alignments = 0
    num_comparisons = 0

    n = len(t)
    m = len(p)

    for i in range(n - m + 1):
        num_alignments += 1
        match = True
        for j in range(m):
            num_comparisons += 1
            if t[i + j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)

    print(f"occurrences: {occurrences}, num_alignments: {num_alignments}, num_comparisons: {num_comparisons}")
    return occurrences, num_alignments, num_comparisons


def boyer_moore_with_counts(p, p_bm, t):
    """
    Boyer–Moore exact matching using preprocessed pattern p_bm (BoyerMoore object).

    Args:
        p: pattern (string)
        p_bm: BoyerMoore(p) object from bm_preproc
        t: text (string)

    Returns:
        occurrences: list of starting indices where p matches t
        num_alignments: how many alignments of p against t were tried
        num_comparisons: how many character comparisons were performed
    """
    occurrences = []
    num_alignments = 0
    num_comparisons = 0

    n = len(t)
    m = len(p)
    i = 0  # current alignment: p[0] aligned with t[i]

    while i <= n - m:
        num_alignments += 1
        shift = 1          # minimum shift if nothing else applies
        mismatched = False

        # compare from right to left
        for j in range(m - 1, -1, -1):
            num_comparisons += 1
            if p[j] != t[i + j]:
                # mismatch at pattern index j and text index i+j
                bad_char_shift = p_bm.bad_character_rule(j, t[i + j])
                good_suffix_shift = p_bm.good_suffix_rule(j)
                shift = max(shift, bad_char_shift, good_suffix_shift)
                mismatched = True
                break

        if not mismatched:
            # full match at position i
            occurrences.append(i)
            # use good-suffix match skip (may be zero)
            match_shift = p_bm.match_skip()
            shift = max(shift, match_shift)

        i += shift
    print(f"occurrences: {occurrences}, num_alignments: {num_alignments}, num_comparisons: {num_comparisons}")
    return occurrences, num_alignments, num_comparisons


def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome





def approximate_match_pigeonhole(P, T, index, max_mismatches=2):
    """
    """
    assert len(P) == 24, "Pattern P must have length 24"
    assert index.k == 8, "Index must be built with k = 8"

    k = index.k
    m = len(P)
    num_parts = m // k  # 3 pieces of length 8

    matches = set()  # avoid duplicates if found via multiple chunks

    # Split P into 3 blocks of 8 bases: P[0:8], P[8:16], P[16:24]
    for part in range(num_parts):
        # Offset of this block inside P
        p_offset = part * k
        P_block = P[p_offset:p_offset + k]

        # Use index to find exact matches of this block
        hits = index.query(P_block)

        for h in hits:
            # h is the position in T where this 8-mer (P_block) begins
            # If P_block starts at h in T, then the whole P would start at:
            start = h - p_offset

            # Check that this proposed alignment is within bounds of T
            if start < 0 or start + m > len(T):
                continue

            # Verify full alignment: count mismatches over all 24 characters
            mismatches = 0
            for j in range(m):
                if T[start + j] != P[j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break

            if mismatches <= max_mismatches:
                matches.add(start)

    return sorted(matches)


def approximate_match_pigeonhole_with_hits(P, T, index, max_mismatches=2):
    total_hits = 0
    matches = set()
    k = index.k
    m = len(P)
    
    for part in range(3):  # 3 blocks of 8
        p_offset = part * k
        block = P[p_offset:p_offset + k]
        
        hits = index.query(block)
        total_hits += len(hits)   # <-- count ALL index hits
        
        for h in hits:
            start = h - p_offset
            if start < 0 or start + m > len(T):
                continue
            mismatches = sum(1 for j in range(m) if P[j] != T[start+j])
            if mismatches <= max_mismatches:
                matches.add(start)

    print(f"length of mathches: {len(matches)}, total hits: {total_hits}")
    return sorted(matches), total_hits



if __name__ == "__main__":
    # tiny demo
    P = "TAATAAA"
    T = "GGCTATAATGCGTA"
    bm = BoyerMoore(P, alphabet='ACGT')
    occ_n, align_n, comp_n = naive_with_counts(P, T)
    occ_bm, align_bm, comp_bm = boyer_moore_with_counts(P, bm, T)

    print("Naive:")
    print("  occurrences     =", occ_n)
    print("  alignments tried =", align_n)
    print("  comparisons      =", comp_n)

    print("\nBoyer–Moore:")
    print("  occurrences     =", occ_bm)
    print("  alignments tried =", align_bm)
    print("  comparisons      =", comp_bm)
