#!/usr/bin/env python3
import argparse
import os
import sys

try:
    import pandas as pd
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: This script requires 'pandas' and 'matplotlib'.")
    print("Please install them using: pip install pandas matplotlib")
    sys.exit(1)


def plot_benchmarks(csv_path, output_dir):
    if not os.path.exists(csv_path):
        print(f"Error: CSV file '{csv_path}' not found.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"[+] Created output directory: {output_dir}")

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    if df.empty:
        print("Error: CSV file is empty.")
        return

    # Clean genome names (remove extension)
    df["Genome"] = df["Genome"].apply(
        lambda x: os.path.splitext(os.path.basename(x))[0]
    )

    # Define genome order by length (approximate)
    # Adjust this list based on your actual files
    genome_order = [
        "phix174",
        "lambda_phage",
        "sars_cov2",
        "yeast_chrI",
        "mycoplasma_genitalium",
        "e_coli_k12",
    ]

    # Filter out genomes not in our known list if necessary, or just sort what we have
    # Create a categorical type for ordering
    df["Genome"] = pd.Categorical(df["Genome"], categories=genome_order, ordered=True)
    df = df.sort_values("Genome")

    # Calculate Search Time per Read (microseconds)
    # This is important because Brute Force uses fewer reads than others
    df["Search Time per Read (us)"] = (df["Search Time (s)"] / df["Num Reads"]) * 1e6

    # Set plot style
    plt.style.use("ggplot")

    # --- Plot 1: Preprocessing Time (Line Chart) ---
    plt.figure(figsize=(10, 6))
    pivot_preproc = df.pivot(
        index="Genome", columns="Algorithm", values="Preproc Time (s)"
    )
    # Reindex to ensure correct order on X-axis
    pivot_preproc = pivot_preproc.reindex(genome_order)

    for column in pivot_preproc.columns:
        plt.plot(pivot_preproc.index, pivot_preproc[column], marker="o", label=column)

    plt.title("Preprocessing Time by Genome Size")
    plt.ylabel("Time (seconds)")
    plt.xlabel("Genome (Increasing Size)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "preprocessing_time_line.png"))
    plt.close()
    print(f"[+] Saved preprocessing_time_line.png")

    # --- Plot 2: Search Time per Read (Line Chart) ---
    plt.figure(figsize=(10, 6))
    pivot_search = df.pivot(
        index="Genome", columns="Algorithm", values="Search Time per Read (us)"
    )
    pivot_search = pivot_search.reindex(genome_order)

    for column in pivot_search.columns:
        plt.plot(pivot_search.index, pivot_search[column], marker="o", label=column)

    plt.title("Search Time per Read (Normalized)")
    plt.ylabel("Time (microseconds)")
    plt.xlabel("Genome (Increasing Size)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "search_time_line.png"))
    plt.close()
    print(f"[+] Saved search_time_line.png")

    # --- Plot 3: Peak Memory Usage (Line Chart) ---
    if df["Peak Memory (MB)"].sum() > 0:
        plt.figure(figsize=(10, 6))
        pivot_memory = df.pivot(
            index="Genome", columns="Algorithm", values="Peak Memory (MB)"
        )
        pivot_memory = pivot_memory.reindex(genome_order)

        for column in pivot_memory.columns:
            plt.plot(pivot_memory.index, pivot_memory[column], marker="o", label=column)

        plt.title("Peak Memory Usage")
        plt.ylabel("Memory (MB)")
        plt.xlabel("Genome (Increasing Size)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "peak_memory_line.png"))
        plt.close()
        print(f"[+] Saved peak_memory_line.png")
    else:
        print("[-] Skipping memory plot (no memory data found).")

    # --- Plot 4: Total Time (Log Scale) ---
    # Useful if Brute Force is included and dominates the scale
    plt.figure(figsize=(12, 6))
    pivot_total = df.pivot(index="Genome", columns="Algorithm", values="Total Time (s)")
    pivot_total.plot(kind="bar", ax=plt.gca())
    plt.title("Total Execution Time (Log Scale)")
    plt.ylabel("Time (seconds) - Log Scale")
    plt.yscale("log")
    plt.xlabel("Genome")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "total_time_log.png"))
    plt.close()
    print(f"[+] Saved total_time_log.png")


def main():
    parser = argparse.ArgumentParser(description="Generate plots from benchmark CSV.")
    parser.add_argument("--csv", required=True, help="Path to benchmark results CSV")
    parser.add_argument("--output", default="./plots", help="Directory to save plots")

    args = parser.parse_args()

    plot_benchmarks(args.csv, args.output)


if __name__ == "__main__":
    main()
