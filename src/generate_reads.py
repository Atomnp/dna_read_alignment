#!/usr/bin/env python3
import argparse
import random
import sys
import os

# Ensure we can import from the same directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from utils import read_fasta


def generate_reads(sequence, num_reads, read_length):
    """
    Generate n random reads of fixed length from the sequence.
    """
    genome_len = len(sequence)
    if genome_len < read_length:
        print(
            f"Error: Genome length ({genome_len}) is smaller than read length ({read_length})"
        )
        sys.exit(1)

    reads = []
    for i in range(num_reads):
        start_pos = random.randint(0, genome_len - read_length)
        read = sequence[start_pos : start_pos + read_length]
        reads.append(read)
    return reads


def save_reads_fasta(reads, output_path):
    """
    Save reads to a file in FASTA format.
    """
    try:
        with open(output_path, "w") as f:
            for i, read in enumerate(reads):
                f.write(f">read_{i}\n")
                f.write(f"{read}\n")
    except Exception as e:
        print(f"Error saving reads: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Generate simulated reads from a reference genome."
    )
    parser.add_argument(
        "--fasta", required=True, help="Path to input reference FASTA file"
    )
    parser.add_argument(
        "--num_reads", type=int, required=True, help="Number of reads to generate"
    )
    parser.add_argument(
        "--read_len", type=int, default=100, help="Length of each read (default: 100)"
    )
    parser.add_argument(
        "--output", required=True, help="Path to output FASTA file for generated reads"
    )

    args = parser.parse_args()

    print(f"[+] Reading reference genome: {args.fasta}")
    seq = read_fasta(args.fasta)
    print(f"[+] Genome length: {len(seq)}")

    print(f"[+] Generating {args.num_reads} reads of length {args.read_len}...")
    reads = generate_reads(seq, args.num_reads, args.read_len)

    print(f"[+] Saving reads to {args.output}...")
    save_reads_fasta(reads, args.output)
    print("[+] Done.")


if __name__ == "__main__":
    main()
