#!/usr/bin/env python3
import argparse
import time
from utils import read_fasta

import argparse
import time
from utils import read_fasta

def brute_force_search(text, pattern):
    """
    Find all occurrences of a pattern in a text using brute-force search.
    """
    n = len(text)
    m = len(pattern)
    matches = []
    for i in range(n - m + 1):
        if text[i:i+m] == pattern:
            matches.append(i)
    return matches

def main():
    parser = argparse.ArgumentParser(
        description="Brute-force Pattern Search in a FASTA File"
    )
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--pattern", required=True, help="Pattern string to search")
    args = parser.parse_args()

    # Time FASTA reading
    print("[+] Loading FASTA sequence...")
    start_time = time.time()
    seq = read_fasta(args.fasta)
    end_time = time.time()
    print(f"[+] Loaded sequence of length: {len(seq)}")
    print(f"[+] Time to load FASTA: {end_time - start_time:.4f} seconds")

    pattern = args.pattern.upper()
    print(f"[+] Searching for: {pattern}")

    # Time the search
    start_time = time.time()
    matches = brute_force_search(seq, pattern)
    end_time = time.time()
    print(f"[+] Time for brute-force search: {end_time - start_time:.4f} seconds")

    if matches:
        print(f"[+] Found {len(matches)} match(es): {matches}")
    else:
        print("[+] No matches found.")


if __name__ == "__main__":
    main()
