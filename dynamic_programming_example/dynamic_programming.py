


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def best_edit_match(P, T):
    m, n = len(P), len(T)
    D = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1):
        D[i][0] = i         # cost to match P[:i] to empty
    for j in range(n+1):
        D[0][j] = 0         # <-- end-free: allow match to start anywhere in T

    for i in range(1, m+1):
        for j in range(1, n+1):
            sub = D[i-1][j-1] + (P[i-1] != T[j-1])
            delete = D[i-1][j] + 1
            insert = D[i][j-1] + 1
            D[i][j] = min(sub, delete, insert)

    # best end position and distance
    j_best = min(range(n+1), key=lambda j: D[m][j])
    best_dist = D[m][j_best]

    # backtrack to find start of the match in T
    i, j = m, j_best
    while i > 0 and j > 0:
        if D[i][j] == D[i-1][j-1] + (P[i-1] != T[j-1]):
            i, j = i-1, j-1
        elif D[i][j] == D[i-1][j] + 1:
            i -= 1
        else:
            j -= 1
    start, end = j, j_best
    return best_dist, start, end, T[start:end]


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


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

from collections import defaultdict

def build_kmer_index(reads, k):
    """Map each k-mer to the set of read indices that CONTAIN it (any position)."""
    idx = defaultdict(set)
    for i, r in enumerate(reads):
        if len(r) < k:
            continue
        for pos in range(len(r) - k + 1):
            idx[r[pos:pos+k]].add(i)
    return idx

def find_overlap_edges(reads, k=30):
    """
    Return a set of directed edges (i, j) where reads[i]'s suffix overlaps
    reads[j]'s prefix by at least k (exact match). Ignores self-overlaps.
    """
    idx = build_kmer_index(reads, k)
    edges = set()

    for i, a in enumerate(reads):
        if len(a) < k:
            continue
        suf = a[-k:]                      # length-k suffix of a
        candidates = idx.get(suf, set())
        for j in candidates:
            if j == i:
                continue                  # no self-overlap
            b = reads[j]
            # verify full suffix/prefix match length >= k
            if overlap(a, b, min_length=k) > 0:
                edges.add((i, j))
    return edges

def count_overlap_edges_fastq(fastq_path, k=30):
    sequences, _ = readFastq(fastq_path)   # names/qualities ignored
    edges = find_overlap_edges(sequences, k)
    return len(edges)

def count_nodes_with_outgoing_fastq(fastq_path, k=30):
    sequences, _ = readFastq(fastq_path)
    edges = find_overlap_edges(sequences, k)
    # nodes that appear as a source at least once
    sources = {i for (i, _) in edges}
    return len(sources)