import sys
from Bio.Blast import NCBIWWW, NCBIXML

# --- Step 1. Get sequence from command line ---
if len(sys.argv) < 2:
    print("Usage: python blast_species_finder.py \"AGCTGACTGCAAGTGACAGTCAGTGAAGT\"")
    sys.exit(1)

sequence = sys.argv[1]

# --- Step 2. Run BLAST search over the internet ---
print("Running BLAST search... please wait (this may take a minute)")
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

# --- Step 3. Parse BLAST results ---
blast_record = NCBIXML.read(result_handle)

# --- Step 4. Display top 3 hits with species information ---
if blast_record.alignments:
    print("\nTop BLAST hits:")
    for i, alignment in enumerate(blast_record.alignments[:3], start=1):
        print(f"{i}. {alignment.title}")
else:
    print("No significant matches found!")
