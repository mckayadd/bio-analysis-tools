#!/usr/bin/env python3
"""
boyer_moore_with_skips.py

Simple Boyer–Moore implementation for teaching:
- Takes T (text) and P (pattern) from the command line
- Locates all occurrences of P in T
- Counts how many "alignments" are skipped thanks to
  - bad character rule
  - good suffix rule

Definition of "skips":
If a rule suggests shifting by S positions, that means it skips (S - 1)
possible alignments of the pattern against the text. We accumulate
these per rule across the full search.
"""

import argparse


# ---------- Bad character preprocessing ----------

def build_bad_char_table(pattern: str):
    """
    Build the "last occurrence" table for the bad character rule.

    For each character c in the pattern, store the index of its
    rightmost occurrence. Characters not in the pattern are treated
    as having index -1 during the search.
    """
    last = {}
    for i, c in enumerate(pattern):
        last[c] = i
    return last


# ---------- Good suffix preprocessing ----------

def build_good_suffix_table(pattern: str):
    """
    Build the good suffix shift table using the standard Boyer–Moore
    preprocessing.

    Returns:
        shift: list of length m+1
        border_pos: list of length m+1 (not needed in search, but
                    useful to keep around for clarity)
    shift[j] is the amount we can shift when a mismatch occurs
    at position j-1 (0-based), i.e. we use shift[j] when mismatch
    is at index (j-1) in the pattern.

    We also use shift[0] when we have found a full match.
    """
    m = len(pattern)
    shift = [0] * (m + 1)
    border_pos = [0] * (m + 1)

    i = m
    j = m + 1
    border_pos[i] = j

    # First pass: compute border positions and some shifts
    while i > 0:
        # Move j back until we find a border that can be extended
        while j <= m and pattern[i - 1] != pattern[j - 1]:
            if shift[j] == 0:
                shift[j] = j - i
            j = border_pos[j]
        i -= 1
        j -= 1
        border_pos[i] = j

    # Second pass: fill in zero shifts with widest border
    j = border_pos[0]
    for i in range(m + 1):
        if shift[i] == 0:
            shift[i] = j
        if i == j:
            j = border_pos[j]

    return shift, border_pos


# ---------- Boyer–Moore search ----------

def boyer_moore_search(text: str, pattern: str):
    """
    Run Boyer–Moore on `text` with `pattern`.

    Returns:
        positions: list of starting indices where pattern is found
        bc_skips: total alignments skipped by bad character rule
        gs_skips: total alignments skipped by good suffix rule
    """
    n = len(text)
    m = len(pattern)
    if m == 0:
        return [], 0, 0, 0

    last = build_bad_char_table(pattern)
    shift, border_pos = build_good_suffix_table(pattern)

    positions = []
    visited_alignments = set()   # which i values we actually check

    bc_skips = 0    # alignments skipped where bad character ruled
    gs_skips = 0    # alignments skipped where good suffix ruled

    total_positions = n - m + 1  # total possible alignments
    i = 0

    while i <= n - m:
        visited_alignments.add(i)

        j = m - 1
        while j >= 0 and pattern[j] == text[i + j]:
            j -= 1

        if j < 0:
            # full match
            positions.append(i)
            chosen_shift = shift[0]  # good suffix shift after a match
            # how many alignments are actually skipped (don't overshoot)
            remaining = max(0, total_positions - (i + 1))
            actual_skip = min(chosen_shift - 1, remaining)
            gs_skips += actual_skip
            i += chosen_shift
        else:
            bad_char = text[i + j]
            last_pos = last.get(bad_char, -1)
            bc_shift = max(1, j - last_pos)
            gs_shift = shift[j + 1]

            chosen_shift = max(bc_shift, gs_shift)
            remaining = max(0, total_positions - (i + 1))
            actual_skip = min(chosen_shift - 1, remaining)

            if bc_shift > gs_shift:
                bc_skips += actual_skip
            elif gs_shift > bc_shift:
                gs_skips += actual_skip
            else:
                # tie – you could attribute to either, or split
                bc_skips += actual_skip

            i += chosen_shift

    total_skipped = total_positions - len(visited_alignments)
    return positions, bc_skips, gs_skips, total_skipped


# ---------- Command-line interface ----------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Boyer–Moore search with skip counting.\n"
            "Example:\n"
            "  python boyer_moore_with_skips.py GGCTATAATGCGTA TAATAAA"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("text", help="T: the text (e.g. genome or long sequence)")
    parser.add_argument("pattern", help="P: the pattern to search for")
    args = parser.parse_args()

    T = args.text
    P = args.pattern

    positions, bc_skips, gs_skips, total_skipped = boyer_moore_search(T, P)

    if positions:
        print(f"Pattern found at positions (0-based): {positions}")
    else:
        print("Pattern not found.")

    print(f"Total bad character skips: {bc_skips}")
    print(f"Total good suffix skips: {gs_skips}")
    print(f"Total skipped: {total_skipped}")


if __name__ == "__main__":
    main()
