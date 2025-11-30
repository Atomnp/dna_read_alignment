#!/usr/bin/env python3
import argparse
import sys
import time
from utils import read_fasta


def get_buckets(s):
    counts = {}
    for char in s:
        counts[char] = counts.get(char, 0) + 1
    return counts


def get_bucket_starts(counts, alphabet):
    starts = {}
    total = 0
    for char in alphabet:
        starts[char] = total
        total += counts[char]
    return starts


def get_bucket_ends(counts, alphabet):
    ends = {}
    total = 0
    for char in alphabet:
        total += counts[char]
        ends[char] = total
    return ends


def sa_is(s):
    """
    Construct suffix array using SA-IS algorithm.
    """
    # Map string to integers for the algorithm
    alphabet = sorted(list(set(s)))
    char_map = {c: i + 1 for i, c in enumerate(alphabet)}
    s_int = [char_map[c] for c in s] + [0]  # Add sentinel

    sa = sa_is_core(s_int, len(alphabet) + 1)
    return sa[1:]  # Remove sentinel index


def sa_is_core(s, K):
    n = len(s)
    t = [False] * n  # False: S-type, True: L-type

    # Classify suffixes
    # t[n-1] is S-type (sentinel is smallest)
    # S-type: s[i] < s[i+1] or (s[i] == s[i+1] and t[i+1] is S)
    # L-type: s[i] > s[i+1] or (s[i] == s[i+1] and t[i+1] is L)
    # Here we use False for S, True for L to match typical implementations or vice versa.
    # Let's stick to: S=False, L=True.
    # Sentinel is S-type.

    # Actually, standard is:
    # S-type: s[i] < s[i+1] or (s[i] == s[i+1] and s[i+1] is S)
    # L-type: s[i] > s[i+1] or (s[i] == s[i+1] and s[i+1] is L)
    # Sentinel is always S-type.

    # Let's use: S=0, L=1.
    # t[i] = 1 if L-type, 0 if S-type
    t = [0] * n
    t[n - 1] = 0  # Sentinel is S-type
    for i in range(n - 2, -1, -1):
        if s[i] < s[i + 1]:
            t[i] = 0  # S
        elif s[i] > s[i + 1]:
            t[i] = 1  # L
        else:
            t[i] = t[i + 1]

    # LMS: Leftmost S-type. t[i] is S and t[i-1] is L.
    is_lms = [False] * n
    lms_indices = []
    for i in range(1, n):
        if t[i] == 0 and t[i - 1] == 1:
            is_lms[i] = True
            lms_indices.append(i)

    # Induced Sort
    sa = [-1] * n
    induced_sort(s, sa, t, is_lms, lms_indices, K)

    # Rename LMS substrings
    new_lms_indices = [i for i in sa if is_lms[i]]

    name = 0
    names = [-1] * n
    names[new_lms_indices[0]] = 0
    prev = new_lms_indices[0]

    for i in range(1, len(new_lms_indices)):
        curr = new_lms_indices[i]
        # Compare LMS substrings
        diff = False
        for j in range(n):
            if s[curr + j] != s[prev + j] or t[curr + j] != t[prev + j]:
                diff = True
                break
            if j > 0 and (is_lms[curr + j] or is_lms[prev + j]):
                break

        if diff:
            name += 1
        names[curr] = name
        prev = curr

    # Recursive step
    s1 = [names[i] for i in lms_indices]
    if name < len(lms_indices) - 1:
        sa1 = sa_is_core(s1, name + 1)
    else:
        # Directly calculate SA if names are unique
        sa1 = [0] * len(s1)
        for i, x in enumerate(s1):
            sa1[x] = i

    # Induce SA from sorted LMS
    sorted_lms = [lms_indices[i] for i in sa1]
    sa = [-1] * n
    induced_sort(s, sa, t, is_lms, sorted_lms, K)

    return sa


def induced_sort(s, sa, t, is_lms, lms_indices, K):
    n = len(s)
    buckets = get_buckets(s)
    alphabet = sorted(buckets.keys())

    # 1. Place LMS characters at end of S buckets
    bucket_ends = get_bucket_ends(buckets, alphabet)
    sa[:] = [-1] * n

    for i in reversed(lms_indices):
        char = s[i]
        pos = bucket_ends[char] - 1
        sa[pos] = i
        bucket_ends[char] -= 1

    # 2. Induce L-types
    bucket_starts = get_bucket_starts(buckets, alphabet)
    for i in range(n):
        if sa[i] > 0:
            j = sa[i] - 1
            if t[j] == 1:  # L-type
                char = s[j]
                pos = bucket_starts[char]
                sa[pos] = j
                bucket_starts[char] += 1

    # 3. Induce S-types
    bucket_ends = get_bucket_ends(buckets, alphabet)
    for i in range(n - 1, -1, -1):
        if sa[i] > 0:
            j = sa[i] - 1
            if t[j] == 0:  # S-type
                char = s[j]
                pos = bucket_ends[char] - 1
                sa[pos] = j
                bucket_ends[char] -= 1


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
        description="SA-IS Suffix Array Construction and Search"
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

    print("[+] Building suffix array using SA-IS...")
    start_time = time.time()
    sa = sa_is(seq)
    end_time = time.time()
    print("[+] Suffix array built.")
    print(
        f"[+] Time to build suffix array (SA-IS): {end_time - start_time:.4f} seconds"
    )

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
