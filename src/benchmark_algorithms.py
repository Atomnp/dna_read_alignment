#!/usr/bin/env python3
import argparse
import time
import sys
import os
import csv
import tracemalloc
from utils import read_fasta

# Ensure we can import from the same directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import brute_force
import suffix_array
import sa_is
import fm_index


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


def benchmark_brute_force(genome, reads, measure_memory=False):
    print("  Running Brute Force...")
    if measure_memory:
        tracemalloc.start()
    start_time = time.perf_counter()

    matches_count = 0
    for read in reads[:200]:
        # print("Searching for read:", read)
        matches = brute_force.brute_force_search(genome, read)
        matches_count += len(matches)

    end_time = time.perf_counter()
    peak = 0
    if measure_memory:
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


def benchmark_naive_sa(genome, reads, measure_memory=False):
    print("  Running Naive Suffix Array...")

    # Preprocessing (Build SA)
    if measure_memory:
        tracemalloc.start()
    start_time = time.perf_counter()
    sa = suffix_array.build_suffix_array(genome)
    total_bytes = deep_getsizeof_list(sa)

    print(f"Suffix Array Memory: {total_bytes / (1024*1024):.2f} MB")
    end_time = time.perf_counter()
    peak_preprocess = 0
    if measure_memory:
        _, peak_preprocess = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    preprocess_time = end_time - start_time

    # Search
    if measure_memory:
        tracemalloc.start()
    start_time = time.perf_counter()
    matches_count = 0
    for read in reads:
        matches = suffix_array.suffix_array_search(genome, sa, read)
        matches_count += len(matches)
    end_time = time.perf_counter()
    peak_search = 0
    if measure_memory:
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


def benchmark_sa_is(genome, reads, sa_save_path=None, measure_memory=False):
    print("  Running SA-IS...")

    # Preprocessing (Build SA)
    if measure_memory:
        tracemalloc.start()
    start_time = time.perf_counter()
    sa = sa_is.sa_is(genome)
    total_bytes = deep_getsizeof_list(sa)
    print(f"Suffix Array Memory: {total_bytes / (1024*1024):.2f} MB")
    end_time = time.perf_counter()
    peak_preprocess = 0
    if measure_memory:
        _, peak_preprocess = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    preprocess_time = end_time - start_time

    # Optionally save SA to file
    if sa_save_path:
        save_suffix_array(sa, sa_save_path)

    # Search
    if measure_memory:
        tracemalloc.start()
    start_time = time.perf_counter()
    matches_count = 0
    for read in reads:
        matches = sa_is.suffix_array_search(genome, sa, read)
        matches_count += len(matches)
    end_time = time.perf_counter()
    peak_search = 0
    if measure_memory:
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


def benchmark_fm_index(genome, reads, measure_memory=False):
    print("  Running FM-Index...")

    # Preprocessing (Build FM-Index)
    if measure_memory:
        tracemalloc.start()
    start_time = time.perf_counter()
    fm = fm_index.FMIndex(genome)

    # Calculate memory of FM Index structure roughly
    total_bytes = deep_getsizeof_list(fm.sa)
    total_bytes += sys.getsizeof(fm.bwt)
    total_bytes += sys.getsizeof(fm.C)
    for char in fm.occ:
        total_bytes += deep_getsizeof_list(fm.occ[char])

    print(f"FM-Index Memory: {total_bytes / (1024*1024):.2f} MB")

    end_time = time.perf_counter()
    peak_preprocess = 0
    if measure_memory:
        _, peak_preprocess = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    preprocess_time = end_time - start_time

    # Search
    if measure_memory:
        tracemalloc.start()
    start_time = time.perf_counter()
    matches_count = 0
    for read in reads:
        matches = fm.locate(read)
        matches_count += len(matches)
    end_time = time.perf_counter()
    peak_search = 0
    if measure_memory:
        _, peak_search = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    search_time = end_time - start_time

    return {
        "name": "FM-Index",
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


def run_benchmark_multiple_genome(
    genome_dir, reads_dir, output_csv, measure_memory=False, skip_brute=False
):
    if not os.path.exists(genome_dir):
        print(f"Error: Genome directory '{genome_dir}' not found.")
        return
    if not os.path.exists(reads_dir):
        print(f"Error: Reads directory '{reads_dir}' not found.")
        return

    # Find genomes
    genomes = [
        f
        for f in os.listdir(genome_dir)
        if f.lower().endswith((".fasta", ".fna", ".fa"))
    ]
    if not genomes:
        print("No genomes found.")
        return

    # Prepare CSV
    fieldnames = [
        "Genome",
        "Algorithm",
        "Preproc Time (s)",
        "Search Time (s)",
        "Total Time (s)",
        "Peak Memory (MB)",
        "Matches",
        "Num Reads",
        "bp",
    ]

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for genome_file in genomes:
            print(f"\n[+] Processing Genome: {genome_file}")
            genome_path = os.path.join(genome_dir, genome_file)
            try:
                genome = read_fasta(genome_path)
            except Exception as e:
                print(f"  Error reading genome: {e}")
                continue

            genome_length = len(genome)
            base_name = os.path.splitext(genome_file)[0]

            # Define read files
            reads_500_file = f"{base_name}_reads_500.fasta"
            reads_200k_file = f"{base_name}_reads_200k.fasta"

            reads_500_path = os.path.join(reads_dir, reads_500_file)
            reads_200k_path = os.path.join(reads_dir, reads_200k_file)

            # 1. Brute Force (500 reads)
            if os.path.exists(reads_500_path) and not skip_brute:
                print(f"  Loading 500 reads from {reads_500_file}")
                reads_500 = load_reads(reads_500_path)
                res = benchmark_brute_force(
                    genome, reads_500, measure_memory=measure_memory
                )
                row = {
                    "Genome": genome_file,
                    "Algorithm": res["name"],
                    "Preproc Time (s)": res["preprocess_time"],
                    "Search Time (s)": res["search_time"],
                    "Total Time (s)": res["total_time"],
                    "Peak Memory (MB)": res["peak_memory_mb"],
                    "Matches": res["matches"],
                    "Num Reads": len(reads_500),
                }
                writer.writerow(row)
                csvfile.flush()  # Ensure write
            else:
                print(f"Skipping Brute Force.")

            # 2. Others (200k reads)
            if os.path.exists(reads_200k_path):
                print(f"  Loading 200k reads from {reads_200k_file}")
                reads_200k = load_reads(reads_200k_path)

                # SA-IS
                res = benchmark_sa_is(genome, reads_200k, measure_memory=measure_memory)
                row = {
                    "Genome": genome_file,
                    "Algorithm": res["name"],
                    "Preproc Time (s)": res["preprocess_time"],
                    "Search Time (s)": res["search_time"],
                    "Total Time (s)": res["total_time"],
                    "Peak Memory (MB)": res["peak_memory_mb"],
                    "Matches": res["matches"],
                    "Num Reads": len(reads_200k),
                    "bp": genome_length,
                }
                writer.writerow(row)

                # Naive SA
                res = benchmark_naive_sa(
                    genome, reads_200k, measure_memory=measure_memory
                )
                row = {
                    "Genome": genome_file,
                    "Algorithm": res["name"],
                    "Preproc Time (s)": res["preprocess_time"],
                    "Search Time (s)": res["search_time"],
                    "Total Time (s)": res["total_time"],
                    "Peak Memory (MB)": res["peak_memory_mb"],
                    "Matches": res["matches"],
                    "Num Reads": len(reads_200k),
                    "bp": genome_length,
                }
                writer.writerow(row)

                # FM-Index
                # res = benchmark_fm_index(
                #     genome, reads_200k, measure_memory=measure_memory
                # )
                # row = {
                #     "Genome": genome_file,
                #     "Algorithm": res["name"],
                #     "Preproc Time (s)": res["preprocess_time"],
                #     "Search Time (s)": res["search_time"],
                #     "Total Time (s)": res["total_time"],
                #     "Peak Memory (MB)": res["peak_memory_mb"],
                #     "Matches": res["matches"],
                #     "Num Reads": len(reads_200k),
                #     "bp": genome_length,
                # }
                # writer.writerow(row)
                csvfile.flush()
            else:
                print(
                    f"  Warning: {reads_200k_file} not found. Skipping SA/FM benchmarks."
                )

    print(f"\n[+] Benchmark complete. Results saved to {output_csv}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark Genome Mapping Algorithms")
    parser.add_argument("--fasta", help="Path to reference genome FASTA (Single Mode)")
    parser.add_argument("--reads", help="Path to reads FASTA (Single Mode)")
    parser.add_argument(
        "--skip-brute",
        action="store_true",
        help="Skip brute force (slow for large genomes)",
    )
    parser.add_argument(
        "--measure-memory",
        action="store_true",
        help="Enable memory benchmarking using tracemalloc (can be slow)",
    )

    # Batch mode arguments
    parser.add_argument(
        "--genome_dir", help="Directory containing reference genomes (Batch Mode)"
    )
    parser.add_argument(
        "--reads_dir", help="Directory containing simulated reads (Batch Mode)"
    )
    parser.add_argument("--output_csv", help="Path to output CSV file (Batch Mode)")

    args = parser.parse_args()

    # Check for Batch Mode
    if args.genome_dir and args.reads_dir and args.output_csv:
        run_benchmark_multiple_genome(
            args.genome_dir,
            args.reads_dir,
            args.output_csv,
            args.measure_memory,
            args.skip_brute,
        )
        return

    # Check for Single Mode
    if args.fasta and args.reads:
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
            results.append(
                benchmark_brute_force(genome, reads, measure_memory=args.measure_memory)
            )
        else:
            print("  Skipping Brute Force...")

        # 2. SA-IS
        results.append(
            benchmark_sa_is(
                genome,
                reads,
                sa_save_path="sa_is_suffix_array.txt",
                measure_memory=args.measure_memory,
            )
        )

        # 3. Naive Suffix Array
        results.append(
            benchmark_naive_sa(genome, reads, measure_memory=args.measure_memory)
        )

        # 4. FM-Index
        results.append(
            benchmark_fm_index(genome, reads, measure_memory=args.measure_memory)
        )

        print_results(results)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
