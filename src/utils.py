import sys

def read_fasta(path):
    """
    Read FASTA file and return a single concatenated sequence.
    Only the first sequence is loaded.
    """
    seq = []
    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                seq.append(line)
    except Exception as e:
        print(f"ERROR reading FASTA: {e}")
        sys.exit(1)

    if not seq:
        print("ERROR: FASTA file contains no sequence data.")
        sys.exit(1)

    return "".join(seq)
