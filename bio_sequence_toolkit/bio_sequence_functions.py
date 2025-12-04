import math

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


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def naive_with_rc(p, t):
    p_ = reverseComplement(p)
    if p_ == p:
        return naive(p, t)
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        matchF = True
        matchR = True 
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                matchF = False
                break
        for k in range(len(p_)):  # loop over characters
            if t[i+k] != p_[k]:  # compare characters
                matchR = False
                break
        if matchF or matchR:
            occurrences.append(i)  # all chars matched; record
    return occurrences


# def naive_2mm(p, t):
#     occurrences = []
#     for i in range(len(t) - len(p) + 1):  # loop over alignments
#         match = True
#         mismatchCount = 0
#         for j in range(len(p)):  # loop over characters
#             if t[i+j] != p[j]:  # compare characters
#                 mismatchCount += 1
#             if mismatchCount > 2:  # compare characters
#                 match = False
#                 break
#         if match:
#             occurrences.append(i)  # all chars matched; record
#     return occurrences

def naive_2mm(p, t):
    m, n = len(p), len(t)
    occ = []
    for i in range(n - m + 1):
        mismatches = 0
        for j, a in enumerate(p):
            if t[i + j] != a:
                mismatches += 1
                if mismatches > 2:
                    break
        else:
            # only runs if the for-loop wasn't broken
            occ.append(i)
    return occ



def bad_cycle_index(filename):
    sequences, qualities = readFastq(filename)
    if not qualities:
        return -1  # no data

    # longest quality string across reads
    max_len = max(len(q) for q in qualities)

    # accumulators per cycle (position)
    sums = [0] * max_len
    counts = [0] * max_len
    # accumulate Phred scores per position
    for qual in qualities:
        for i, ch in enumerate(qual):
            q = Phred33toQ(ch)          # ord(ch) - 33
            sums[i] += q
            counts[i] += 1

    # mean quality per position (handle positions with zero coverage)
    means = [
        (s / c) if c > 0 else float("inf")
        for s, c in zip(sums, counts)
    ]

    # index (0-based) of the lowest-mean (worst) cycle
    lowest_index = min(range(len(means)), key=lambda i: means[i])
    return lowest_index


def phred_score(char):
    """Convert an ASCII character (Phred+33 encoding) to Q-score and error probability."""
    q = ord(char) - 33                       # Convert ASCII to Phred Q-score
    p_error = 10 ** (-q / 10)                # Compute probability of error
    return q, p_error


def Phred33toQ(char):
    """Convert ASCII (Phred+33) character to numeric Q score."""
    return ord(char) - 33


def QtoPhred33(Q):
    """Convert numeric Q score to ASCII (Phred+33) character."""
    return chr(Q + 33)