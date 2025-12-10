
import itertools

def overlap(a, b, min_length=1):
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def scs_len(ss):
    """Return the length of a shortest common superstring of ss."""
    best = None
    for perm in itertools.permutations(ss):
        sup_len = len(perm[0])
        for i in range(len(ss) - 1):
            olen = overlap(perm[i], perm[i+1], min_length=1)
            sup_len += len(perm[i+1]) - olen
        if best is None or sup_len < best:
            best = sup_len
    return best



def scs_count(ss):
    """Return (#distinct SCS strings, SCS length, set of those SCS strings)."""
    best_len = None
    best_set = set()
    for perm in itertools.permutations(ss):
        sup = perm[0]
        for i in range(len(ss)-1):
            olen = overlap(perm[i], perm[i+1], min_length=1)
            sup += perm[i+1][olen:]
        L = len(sup)
        if best_len is None or L < best_len:
            best_len = L
            best_set = {sup}
        elif L == best_len:
            best_set.add(sup)
    return len(best_set), best_len, best_set


strings = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]
count, L, scs_set = scs_count(strings)
print("SCS length:", L)              # 11
print("Number of distinct SCS:", count)  # 4
print("Examples:", sorted(scs_set))