#!/usr/bin/env python3
import argparse
import time
from bisect import bisect_left, bisect_right
from utils import read_fasta


def build_suffix_array(s):
    """
    Build suffix array using prefix-doubling method.
    Returns SA: list of starting indices of sorted suffixes.
    """
    n = len(s)
    k = 1
    # rank[i] = rank of suffix at position i
    rank = list(map(ord, s))
    tmp = [0] * n
    sa = list(range(n))

    while True:
        # Sort by pair (rank[i], rank[i+k])
        sa.sort(key=lambda i: (rank[i], rank[i + k] if i + k < n else -1))

        # Recompute temp ranks
        tmp[sa[0]] = 0
        for i in range(1, n):
            prev = sa[i - 1]
            curr = sa[i]
            tmp[curr] = tmp[prev]
            if (rank[curr], rank[curr + k] if curr + k < n else -1) != (
                rank[prev],
                rank[prev + k] if prev + k < n else -1,
            ):
                tmp[curr] += 1

        rank = tmp[:]
        k <<= 1  # double k

        if max(rank) == n - 1:
            break

    return sa


def compare_substring(s, pattern, start):
    """
    Compare pattern with s[start: start+len(pattern)].
    Return:
      -1 if s_sub < pattern
       0 if equal (prefix match)
       1 if s_sub > pattern
    """
    m = len(pattern)
    sub = s[start : start + m]
    if sub == pattern:
        return 0
    return -1 if sub < pattern else 1


def suffix_array_search(s, sa, pattern):
    """
    Return all positions where pattern appears in s.
    Uses two binary searches: lower bound and upper bound.
    """

    # --- lower bound: first index with suffix >= pattern ---
    left, right = 0, len(sa)
    while left < right:
        mid = (left + right) // 2
        if compare_substring(s, pattern, sa[mid]) >= 0:
            right = mid
        else:
            left = mid + 1
    lb = left

    left, right = lb, len(sa)
    while left < right:
        mid = (left + right) // 2
        if compare_substring(s, pattern, sa[mid]) > 0:
            right = mid
        else:
            left = mid + 1
    ub = left

    return sa[lb:ub]


def main():
    parser = argparse.ArgumentParser(
        description="Suffix Array Pattern Search in a FASTA File"
    )
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--pattern", required=True, help="Pattern string to search")
    args = parser.parse_args()

    print("[+] Loading FASTA sequence...")
    start_time = time.time()
    seq = read_fasta(args.fasta)
    end_time = time.time()
    print(f"[+] Loaded sequence of length: {len(seq)}")
    print(f"[+] Time to load FASTA: {end_time - start_time:.4f} seconds")

    print("[+] Building suffix array (this may take time for large genomes)...")
    start_time = time.time()
    sa = build_suffix_array(seq)
    end_time = time.time()
    print("[+] Suffix array built.")
    print(f"[+] Time to build suffix array: {end_time - start_time:.4f} seconds")

    pattern = args.pattern.upper()
    print(f"[+] Searching for: {pattern}")

    start_time = time.time()
    matches = suffix_array_search(seq, sa, pattern)
    end_time = time.time()
    print(f"[+] Time for suffix array search: {end_time - start_time:.4f} seconds")

    if matches:
        print(f"[+] Found {len(matches)} match(es): {matches}")
    else:
        print("[+] No matches found.")


if __name__ == "__main__":
    main()
