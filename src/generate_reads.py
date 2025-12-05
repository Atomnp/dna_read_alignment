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
        description="Generate simulated reads from reference genomes in a directory."
    )
    parser.add_argument(
        "--input_dir",
        default="./genome_data",
        help="Directory containing reference FASTA files (default: ./genome_data)",
    )
    parser.add_argument(
        "--output_dir",
        default="./simulated_reads",
        help="Directory to save simulated reads (default: ./simulated_reads)",
    )
    parser.add_argument(
        "--read_len", type=int, default=100, help="Length of each read (default: 100)"
    )

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"[+] Created output directory: {args.output_dir}")

    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist.")
        sys.exit(1)

    # Find all FASTA files
    files = [
        f
        for f in os.listdir(args.input_dir)
        if f.lower().endswith((".fasta", ".fna", ".fa"))
    ]

    if not files:
        print(f"No FASTA files found in {args.input_dir}")
        sys.exit(0)

    print(f"[+] Found {len(files)} FASTA files in {args.input_dir}")

    for filename in files:
        file_path = os.path.join(args.input_dir, filename)
        print(f"\n[+] Processing: {filename}")

        # read_fasta might sys.exit(1) on error, which stops the script.
        # Assuming valid inputs for now as per utils.py design.
        seq = read_fasta(file_path)
        print(f"    Genome length: {len(seq)}")

        # Configurations: 500 reads and 200k reads
        configs = [(500, "500"), (200000, "200k")]

        base_name = os.path.splitext(filename)[0]

        for num_reads, suffix in configs:
            print(f"    Generating {num_reads} reads...")
            try:
                reads = generate_reads(seq, num_reads, args.read_len)

                output_filename = f"{base_name}_reads_{suffix}.fasta"
                output_path = os.path.join(args.output_dir, output_filename)

                save_reads_fasta(reads, output_path)
                print(f"    Saved: {output_filename}")
            except Exception as e:
                print(f"    Error generating/saving reads: {e}")

    print("\n[+] Batch processing complete.")


if __name__ == "__main__":
    main()
