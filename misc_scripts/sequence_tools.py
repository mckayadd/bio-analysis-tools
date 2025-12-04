import sys
from Bio.Seq import Seq

# --- Step 1. Check for command line argument ---
if len(sys.argv) < 2:
    print("Usage: python sequence_tools.py \"ATGCGTACGTA\"")
    sys.exit(1)

# --- Step 2. Clean and prepare the sequence ---
sequence_str = sys.argv[1].replace("\n", "").replace("\r", "").replace(" ", "").upper()

# --- Step 3. Create Seq object ---
my_seq = Seq(sequence_str)

# --- Step 4. Display DNA info ---
print("Original DNA sequence:       ", my_seq)
print("Reverse complement sequence: ", my_seq.reverse_complement())

# --- Step 5. Translate safely (ignore incomplete codons) ---
try:
    protein = my_seq.translate(to_stop=False)
    print("Translated protein sequence: ", protein)
except Exception as e:
    print("Translation error:", e)
    print("Tip: Make sure the sequence length is a multiple of 3.")


'''
Example run

 % python3 sequence_tools.py "TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTAC
AATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCAC
CTACGGTAGAG"
Original DNA sequence:        TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG
Reverse complement sequence:  CTCTACCGTAGGTGCCAAAGAAGAATGAGAAGAAGCCCGAAAAAGAGCAAAAAAGAACATAACCTCCGATACGATAAACAGAATAAAACCATATCGAGGTCCTAATTGTACGACTTTGGTATGATGTCCTTCGAACGTGGATTCACGTAGAACATCGCGCCACCATACGAACATGGTATATAGGATAAATATGAGGCCCA
/opt/anaconda3/lib/python3.12/site-packages/Bio/Seq.py:2879: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
  warnings.warn(
Translated protein sequence:  WASYLSYIPCSYGGAMFYVNPRSKDIIPKSYN*DLDMVLFCLSYRRLCSFLLFFGLLLILLWHLR*
'''