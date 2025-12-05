#!/usr/bin/env python3
import argparse
import os
import sys

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
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

    # Ensure 'bp' column exists
    if "bp" not in df.columns:
        print("Error: 'bp' column missing in CSV. Cannot plot by genome length.")
        return

    # Sort by bp for consistent plotting
    df = df.sort_values("bp")

    # Set plot style for academic look
    plt.style.use("seaborn-v0_8-whitegrid")

    # Get unique algorithms
    algorithms = df["Algorithm"].unique()

    # --- Plot 1: Indexing Time vs Genome Length (Log-Log) ---
    plt.figure(figsize=(10, 6))

    has_indexing_data = False
    for algo in algorithms:
        subset = df[df["Algorithm"] == algo]
        # Filter out 0 values for log plot (e.g., Brute Force)
        subset = subset[subset["Preproc Time (s)"] > 0]

        if not subset.empty:
            has_indexing_data = True
            plt.plot(
                subset["bp"],
                subset["Preproc Time (s)"],
                marker="o",
                label=algo,
                linewidth=2,
            )

    if has_indexing_data:
        # Annotate genomes (place label at the top-most point for each genome)
        valid_indexing = df[df["Preproc Time (s)"] > 0]
        for genome in valid_indexing["Genome"].unique():
            genome_data = valid_indexing[valid_indexing["Genome"] == genome]
            if not genome_data.empty:
                # Find the point with max Y value to place label above it
                max_idx = genome_data["Preproc Time (s)"].idxmax()
                max_row = genome_data.loc[max_idx]
                plt.annotate(
                    genome,
                    (max_row["bp"], max_row["Preproc Time (s)"]),
                    textcoords="offset points",
                    xytext=(0, 5),
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    fontweight="bold",
                )
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Genome Length (bp)", fontsize=12)
        plt.ylabel("Indexing Time (s)", fontsize=12)
        plt.title("Indexing Time Complexity vs Genome Length", fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(True, which="both", ls="-", alpha=0.3)

        # Format axes
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x:g}"))

        plt.tight_layout()
        save_path = os.path.join(output_dir, "indexing_time_complexity.png")
        plt.savefig(save_path, dpi=300)
        print(f"[+] Saved {save_path}")
        plt.close()
    else:
        print("[-] No non-zero indexing time data found. Skipping plot.")

    # --- Plot 2: Search Time per 100k Queries vs Genome Length (Log-Log) ---
    # Calculate normalized search time
    df["Search Time per 100k"] = (df["Search Time (s)"] / df["Num Reads"]) * 100000

    plt.figure(figsize=(10, 6))

    for algo in algorithms:
        subset = df[df["Algorithm"] == algo]
        plt.plot(
            subset["bp"],
            subset["Search Time per 100k"],
            marker="o",
            label=algo,
            linewidth=2,
        )

    # Annotate genomes
    for genome in df["Genome"].unique():
        genome_data = df[df["Genome"] == genome]
        if not genome_data.empty:
            max_idx = genome_data["Search Time per 100k"].idxmax()
            max_row = genome_data.loc[max_idx]
            plt.annotate(
                genome,
                (max_row["bp"], max_row["Search Time per 100k"]),
                textcoords="offset points",
                xytext=(0, 5),
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Genome Length (bp)", fontsize=12)
    plt.ylabel("Search Time per 100k Queries (s)", fontsize=12)
    plt.title("Search Time Complexity vs Genome Length", fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", ls="-", alpha=0.3)

    # Format axes
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x:g}"))

    plt.tight_layout()
    save_path = os.path.join(output_dir, "search_time_complexity.png")
    plt.savefig(save_path, dpi=300)
    print(f"[+] Saved {save_path}")
    plt.close()

    # --- Plot 3: Total Time (Bar Chart) ---
    # Helper to format bp for labels
    def format_bp(size):
        if size >= 1e9:
            return f"{size/1e9:.1f}B"
        elif size >= 1e6:
            return f"{size/1e6:.1f}M"
        elif size >= 1e3:
            return f"{size/1e3:.1f}K"
        else:
            return str(size)

    df["Size Label"] = df["bp"].apply(format_bp)
    df["X_Label"] = df["Genome"] + "\n(" + df["Size Label"] + ")"

    # Pivot for bar chart
    pivot_total = df.pivot(
        index="X_Label", columns="Algorithm", values="Total Time (s)"
    )
    # Sort index by bp (using the original df to get order)
    # We need to get unique X_Labels in the correct order
    unique_labels = df.sort_values("bp")["X_Label"].unique()
    pivot_total = pivot_total.reindex(unique_labels)

    plt.figure(figsize=(12, 6))
    pivot_total.plot(kind="bar", ax=plt.gca())
    plt.title("Total Execution Time (Log Scale)", fontsize=14)
    plt.ylabel("Time (seconds) - Log Scale", fontsize=12)
    plt.yscale("log")
    plt.xlabel("Genome", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", which="major", alpha=0.3)
    plt.tight_layout()

    save_path = os.path.join(output_dir, "total_time_bar.png")
    plt.savefig(save_path, dpi=300)
    print(f"[+] Saved {save_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Generate academic plots from benchmark CSV."
    )
    parser.add_argument("--csv", required=True, help="Path to benchmark results CSV")
    parser.add_argument("--output", default="./plots", help="Directory to save plots")

    args = parser.parse_args()

    plot_benchmarks(args.csv, args.output)


if __name__ == "__main__":
    main()
