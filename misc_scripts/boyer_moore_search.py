import sys

# --- Helper Functions (Same as before) ---
def build_bad_char_table(pattern):
    m = len(pattern)
    table = {}
    for i in range(m - 1):
        table[pattern[i]] = i
    return table

def build_good_suffix_tables(pattern):
    m = len(pattern)
    border = [0] * m
    j = 0
    for i in range(1, m):
        while j > 0 and pattern[i] != pattern[j]:
            j = border[j - 1]
        if pattern[i] == pattern[j]:
            j += 1
        border[i] = j

    good_suffix_shift = [0] * (m + 1)
    j = 0
    for i in range(m - 1, -1, -1):
        if border[i] == i + 1:
            for k in range(j, m - 1 - i):
                if good_suffix_shift[k] == 0:
                    good_suffix_shift[k] = m - 1 - i
            j = m - 1 - i

    for i in range(m - 1):
        shift = m - 1 - border[i]
        good_suffix_shift[border[i]] = shift

    return good_suffix_shift, m - border[m - 1]

# --- Main Search Function ---
def boyer_moore_total_rule_skips(text, pattern):
    """
    Performs the Boyer-Moore search and aggregates total skips for BC and GS rules.

    Returns:
        tuple: (list of found indices, total BC skips, total GS skips)
    """
    n = len(text)
    m = len(pattern)

    # Pre-processing
    bad_char_table = build_bad_char_table(pattern)
    good_suffix_shifts, matched_prefix_shift = build_good_suffix_tables(pattern)

    i = 0  # current alignment index in T
    found_indices = []
    total_bc_skips = 0
    total_gs_skips = 0

    if m == 0 or n < m:
        return [], 0, 0

    # Searching phase
    while i <= n - m:
        j = m - 1  # index in P, starting from the end
        match_success = True

        # 1. Mismatch detection (right-to-left comparison)
        while j >= 0:
            if pattern[j] != text[i + j]:
                match_success = False

                # --- Bad Character Rule (BC) ---
                char_mismatch = text[i + j]
                last_occurrence_index = bad_char_table.get(char_mismatch, -1)
                bc_shift = j - last_occurrence_index

                # --- Good Suffix Rule (GS) ---
                good_suffix_len = m - 1 - j
                gs_shift = good_suffix_shifts[good_suffix_len] if good_suffix_len > 0 else 1

                # --- Final Shift: Take the maximum and attribute the skip ---
                final_shift = max(bc_shift, gs_shift)
                if final_shift == 0:
                    final_shift = 1 

                if bc_shift >= gs_shift:
                    # BC rule dictated the shift
                    total_bc_skips += final_shift
                else:
                    # GS rule dictated the shift (strictly greater)
                    total_gs_skips += final_shift

                # Perform the skip
                i += final_shift
                break  # Move to the next alignment (outer while loop)

            j -= 1

        # 2. Match success check (j < 0 means the whole pattern matched)
        if match_success:
            found_indices.append(i)

            # Shift after a match (always attributed to the GS/MP rule)
            final_shift = matched_prefix_shift
            total_gs_skips += final_shift # Full match shift is part of the GS/MP logic
            
            i += final_shift

    return found_indices, total_bc_skips, total_gs_skips

if __name__ == '__main__':
    # --- Input Validation ---
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <text_T> <pattern_P>")
        print("Example: python {sys.argv[0]} ABC ABCDABD ABDA")
        sys.exit(1)

    # Get arguments from the command line
    text = sys.argv[1]
    pattern = sys.argv[2]
    m = len(pattern)

    print(f" Text (T):   {text}")
    print(f" Pattern (P): {pattern}")
    print("-" * 30)
    
    # Perform the search and get results
    indices, bc_skips, gs_skips = boyer_moore_total_rule_skips(text, pattern)

    # --- Output Results ---
    print("## Rule-Based Skip Totals")
    print("-" * 30)

    print(f"Bad Character Rule Total Skips: **{bc_skips}**")
    print(f"Good Suffix Rule Total Skips: **{gs_skips}**")
    
    total_skips = bc_skips + gs_skips
    print(f"Combined Total Skips: **{total_skips}**")

    print("\n## Final Results")
    
    if indices:
        print("\n Pattern found at the following starting indices (0-based):")
        for index in indices:
            print(f"   -> Index {index}: {text[:index]}**{text[index:index+m]}**{text[index+m:]}")
    else:
        print("\n Pattern not found in the text.")