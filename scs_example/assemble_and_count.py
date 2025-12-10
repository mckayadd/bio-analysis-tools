import sys
import heapq
import time

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
    """Return length of longest suffix of 'a' matching prefix of 'b'."""
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1


def greedy_scs(reads, k=10, progress_every=50):
    """
    Greedy shortest common superstring using a max-heap of overlaps.
    Uses unique IDs to avoid stale references and index errors.
    """
    # Assign unique IDs to each read
    id_to_seq = {i: s for i, s in enumerate(reads)}
    alive = set(id_to_seq.keys())
    total = len(alive)

    # Build initial max-heap of overlaps: (-overlap, id_a, id_b)
    heap = []
    for i in alive:
        si = id_to_seq[i]
        for j in alive:
            if i == j:
                continue
            sj = id_to_seq[j]
            olen = overlap(si, sj, min_length=k)
            if olen > 0:
                heapq.heappush(heap, (-olen, i, j))

    step = 0
    last_progress_time = time.time()

    while len(alive) > 1:
        step += 1

        # Pop best available overlap that is still valid
        valid_found = False
        while heap:
            neg_olen, i, j = heapq.heappop(heap)
            if i in alive and j in alive:
                olen = -neg_olen
                valid_found = True
                break

        if not valid_found:
            # No overlaps left; concatenate remaining reads to finish
            # This keeps functionality (single assembled string) and avoids stalling.
            remaining = [id_to_seq[rid] for rid in alive]
            assembled = "".join(remaining)
            # Mark all merged to a single ID
            only = min(alive)
            id_to_seq[only] = assembled
            alive = {only}
            # Progress output
            done = total - len(alive)
            percent = (done / total) * 100
            sys.stdout.write(f"\rStep {step}: {len(alive)} reads remaining ({percent:.2f}% complete)")
            sys.stdout.flush()
            break

        # Merge j into i using the overlap
        seq_i = id_to_seq[i]
        seq_j = id_to_seq[j]
        combined = seq_i + seq_j[olen:]
        id_to_seq[i] = combined
        alive.remove(j)  # j is now merged into i

        # Update overlaps only for pairs involving i (the changed sequence)
        for m in list(alive):
            if m == i:
                continue
            sm = id_to_seq[m]
            olen_im = overlap(id_to_seq[i], sm, min_length=k)
            if olen_im > 0:
                heapq.heappush(heap, (-olen_im, i, m))
            olen_mi = overlap(sm, id_to_seq[i], min_length=k)
            if olen_mi > 0:
                heapq.heappush(heap, (-olen_mi, m, i))

        # Progress output (rate-limited)
        if step % progress_every == 0 or (time.time() - last_progress_time) > 1.0:
            last_progress_time = time.time()
            done = total - len(alive)
            percent = (done / total) * 100
            sys.stdout.write(f"\rStep {step}: {len(alive)} reads remaining ({percent:.2f}% complete)")
            sys.stdout.flush()

    print("\nAssembly finished.")
    # Return the single remaining sequence
    final_id = next(iter(alive))
    return id_to_seq[final_id]


def count_AT(genome):
    """Return counts of A and T in genome."""
    return genome.count("A"), genome.count("T")


def main():
    # Configure these if you like, or replace with argparse for CLI usage
    filename = "ads1_week4_reads.fq"  # <-- set your input FASTQ file here
    min_overlap_k = 10                # increase for speed if reads are long
    progress_every = 50               # print progress every N merges

    sequences, qualities = readFastq(filename)
    print(f"Loaded {len(sequences)} reads from {filename}")

    assembled_genome = greedy_scs(sequences, k=min_overlap_k, progress_every=progress_every)
    count_A, count_T = count_AT(assembled_genome)

    print("Assembled genome length:", len(assembled_genome))
    print("Number of A bases:", count_A)
    print("Number of T bases:", count_T)


if __name__ == "__main__":
    main()
