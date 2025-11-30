#!/usr/bin/env python3
import argparse
import time
import sys
import os
import tracemalloc
from utils import read_fasta

# Ensure we can import from the same directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import brute_force
import suffix_array
import sa_is


def load_reads(path):
    reads = []
    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                reads.append(line)
    except Exception as e:
        print(f"Error reading reads file: {e}")
        sys.exit(1)
    return reads


def benchmark_brute_force(genome, reads):
    print("  Running Brute Force...")
    tracemalloc.start()
    start_time = time.perf_counter()

    matches_count = 0
    for read in reads:
        print("Searching for read:", read)
        matches = brute_force.brute_force_search(genome, read)
        matches_count += len(matches)

    end_time = time.perf_counter()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "name": "Brute Force",
        "preprocess_time": 0.0,
        "search_time": end_time - start_time,
        "total_time": end_time - start_time,
        "peak_memory_mb": peak / (1024 * 1024),
        "matches": matches_count,
    }


def deep_getsizeof_list(lst):
    size = sys.getsizeof(lst)  # list object itself
    for item in lst:
        size += sys.getsizeof(item)  # each Python int
    return size


def benchmark_naive_sa(genome, reads):
    print("  Running Naive Suffix Array...")

    # Preprocessing (Build SA)
    tracemalloc.start()
    start_time = time.perf_counter()
    sa = suffix_array.build_suffix_array(genome)
    total_bytes = deep_getsizeof_list(sa)

    print(f"Suffix Array Memory: {total_bytes / (1024*1024):.2f} MB")
    end_time = time.perf_counter()
    _, peak_preprocess = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    preprocess_time = end_time - start_time

    # Search
    tracemalloc.start()
    start_time = time.perf_counter()
    matches_count = 0
    for read in reads:
        matches = suffix_array.suffix_array_search(genome, sa, read)
        matches_count += len(matches)
    end_time = time.perf_counter()
    _, peak_search = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    search_time = end_time - start_time

    return {
        "name": "Naive SA",
        "preprocess_time": preprocess_time,
        "search_time": search_time,
        "total_time": preprocess_time + search_time,
        "peak_memory_mb": max(peak_preprocess, peak_search) / (1024 * 1024),
        "matches": matches_count,
    }


def save_suffix_array(sa, filepath):
    """Save suffix array to a file (one integer per line)."""
    try:
        with open(filepath, "w") as f:
            for idx in sa:
                f.write(f"{idx}\n")
    except Exception as e:
        print(f"Error saving suffix array: {e}")


def load_suffix_array(filepath):
    """Load suffix array from a file (one integer per line)."""
    sa = []
    try:
        with open(filepath, "r") as f:
            for line in f:
                sa.append(int(line.strip()))
    except Exception as e:
        print(f"Error loading suffix array: {e}")
        return None
    return sa


def benchmark_sa_is(genome, reads, sa_save_path=None):
    print("  Running SA-IS...")

    # Preprocessing (Build SA)
    tracemalloc.start()
    start_time = time.perf_counter()
    sa = sa_is.sa_is(genome)
    total_bytes = deep_getsizeof_list(sa)
    print(f"Suffix Array Memory: {total_bytes / (1024*1024):.2f} MB")
    end_time = time.perf_counter()
    _, peak_preprocess = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    preprocess_time = end_time - start_time

    # Optionally save SA to file
    if sa_save_path:
        save_suffix_array(sa, sa_save_path)

    # Search
    tracemalloc.start()
    start_time = time.perf_counter()
    matches_count = 0
    for read in reads:
        matches = sa_is.suffix_array_search(genome, sa, read)
        matches_count += len(matches)
    end_time = time.perf_counter()
    _, peak_search = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    search_time = end_time - start_time

    return {
        "name": "SA-IS",
        "preprocess_time": preprocess_time,
        "search_time": search_time,
        "total_time": preprocess_time + search_time,
        "peak_memory_mb": max(peak_preprocess, peak_search) / (1024 * 1024),
        "matches": matches_count,
    }


def print_results(results):
    print("\n" + "=" * 85)
    print(
        f"{'Algorithm':<15} | {'Preproc (s)':<12} | {'Search (s)':<12} | {'Total (s)':<12} | {'Peak Mem (MB)':<14} | {'Matches':<8}"
    )
    print("-" * 85)
    for res in results:
        print(
            f"{res['name']:<15} | {res['preprocess_time']:<12.4f} | {res['search_time']:<12.4f} | {res['total_time']:<12.4f} | {res['peak_memory_mb']:<14.2f} | {res['matches']:<8}"
        )
    print("=" * 85 + "\n")


def main():
    parser = argparse.ArgumentParser(description="Benchmark Genome Mapping Algorithms")
    parser.add_argument("--fasta", required=True, help="Path to reference genome FASTA")
    parser.add_argument("--reads", required=True, help="Path to reads FASTA")
    parser.add_argument(
        "--skip-brute",
        action="store_true",
        help="Skip brute force (slow for large genomes)",
    )

    args = parser.parse_args()

    print(f"[+] Loading Genome: {args.fasta}")
    genome = read_fasta(args.fasta)
    print(f"    Length: {len(genome)}")
    # print memory usage of genome
    print(f"    Approx. Memory: {sys.getsizeof(genome) / (1024 * 1024):.2f} MB")

    print(f"[+] Loading Reads: {args.reads}")
    reads = load_reads(args.reads)
    print(f"    Count: {len(reads)}")
    # print memory usage of reads
    reads_mem = sum(sys.getsizeof(read) for read in reads)
    print(f"    Approx. Memory: {reads_mem / (1024 * 1024):.2f} MB")

    results = []

    # 1. Brute Force
    if not args.skip_brute:
        results.append(benchmark_brute_force(genome, reads))
    else:
        print("  Skipping Brute Force...")

    # 2. SA-IS
    results.append(
        benchmark_sa_is(genome, reads, sa_save_path="sa_is_suffix_array.txt")
    )

    # 3. Naive Suffix Array
    results.append(benchmark_naive_sa(genome, reads))

    print_results(results)


if __name__ == "__main__":
    main()
