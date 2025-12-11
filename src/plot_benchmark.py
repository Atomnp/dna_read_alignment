#!/usr/bin/env python3
import argparse
import os
import sys

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import numpy as np
except ImportError:
    print("Error: This script requires 'pandas', 'matplotlib', and 'numpy'.")
    print("Please install them using: pip install pandas matplotlib numpy")
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

    # Get unique genomes
    unique_genomes = df["Genome"].unique()

    # --- Plot 1: Indexing Time vs Genome Length (Log-Log) ---
    plt.figure(figsize=(10, 6))

    has_indexing_data = False
    max_time_val = 0

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
            max_time_val = max(max_time_val, subset["Preproc Time (s)"].max())

    min_bp = df["bp"].min()
    max_bp = df["bp"].max()
    x_theory = np.geomspace(min_bp, max_bp, 100)
    y_linear = x_theory
    plt.plot(
        x_theory,
        y_linear,
        color="gray",
        linestyle="-.",
        label="y = x",
        alpha=0.7,
    )

    y_nsq_log_n = x_theory**2 * np.log2(x_theory)
    plt.plot(
        x_theory,
        y_nsq_log_n,
        color="black",
        linestyle="--",
        label="y = n² log(n)",
        alpha=0.7,
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
        algo_name = algo
        if algo == "Brute Force":
            # Skip Brute Force for search time plot due to poor scaling
            algo_name = "Brute Force (Linear scan)"
        elif algo == "Naive SA":
            algo_name = "Naive SA (Binary Search)"
        elif algo == "SA-IS":
            algo_name = "SA-IS (Binary Search)"
        subset = df[df["Algorithm"] == algo]
        plt.plot(
            subset["bp"],
            subset["Search Time per 100k"],
            marker="o",
            label=algo_name,
            linewidth=2,
        )

    # Add theoretical complexity line: y =log2(x)
    min_bp = df["bp"].min()
    max_bp = df["bp"].max()
    x_theory = np.geomspace(min_bp, max_bp, 100)
    y_theory = np.log2(x_theory)
    plt.plot(
        x_theory,
        y_theory,
        color="black",
        linestyle="--",
        label="y = log(x)",
        alpha=0.7,
    )

    y_linear = x_theory
    plt.plot(
        x_theory,
        y_linear,
        color="gray",
        linestyle="-.",
        label="y = x",
        alpha=0.7,
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

    # Plot 4: Final Index Memory vs Genome Size
    fig, ax = plt.subplots(figsize=(10, 6))
    points_plotted = 0

    for algo in algorithms:
        if "Brute" in algo:
            continue

        subset = df[df["Algorithm"] == algo]
        subset = subset[subset["Final Index Memory (MB)"] > 0]

        if not subset.empty:
            ax.plot(
                subset["bp"],
                subset["Final Index Memory (MB)"],
                marker="s",
                label=algo,
                linewidth=2,
            )
            points_plotted += len(subset)

    if points_plotted > 0:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Genome Length (bp)")
        ax.set_ylabel("Final Index Memory (MB)")
        ax.set_title("Memory Scaling: Index Size vs Genome Length")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", ls="-", alpha=0.3)

        # Force tick formatting to avoid scientific notation confusion
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x:g}"))

        plt.tight_layout()
        save_path = os.path.join(output_dir, "memory_scaling.png")
        plt.savefig(save_path, dpi=300)
        print(f"[+] Saved {save_path}")
    else:
        print("[-] Warning: No valid memory data found. Plot skipped.")
    plt.close()

    # Plot 5: Space-Time Trade-off
    num_genomes = len(unique_genomes)
    cols = 2
    rows = math.ceil(num_genomes / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(12, 5 * rows))

    if num_genomes == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for i, genome in enumerate(unique_genomes):
        ax = axes[i]
        subset = df[df["Genome"] == genome]

        for _, row in subset.iterrows():
            algo_name = row["Algorithm"]

            # --- CHANGE: Strict Filter for FM-Index Only ---
            if "FM-Index" not in algo_name:
                continue
            # -----------------------------------------------

            x_val = row["Final Index Memory (MB)"]
            if x_val == 0:
                continue

            y_val = row["Search Time per 100k"]

            # Since we only plot FM-Index, we can keep the blue style or vary it
            color = "tab:blue"
            marker = "o"

            ax.scatter(
                x_val, y_val, s=100, alpha=0.8, c=color, marker=marker, edgecolors="k"
            )

            # Annotate short names
            short = (
                algo_name.split("(")[-1]
                .replace("Checkpoint Rate =", "CP")
                .replace("SA Sampling Rate =", "SA")
                .replace(")", "")
            )
            ax.annotate(
                short,
                (x_val, y_val),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )

        ax.set_title(f"Genome: {genome}", fontsize=11, fontweight="bold")
        ax.set_xlabel("Memory (MB)")
        ax.set_ylabel("Time (s) / 100k Reads")
        ax.grid(True, linestyle="--", alpha=0.5)

        # Use log scale for Y axis (Time)
        ax.set_yscale("log")

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Simplified Legend for just FM-Index
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="tab:blue",
            label="FM-Index Configuration",
            markersize=10,
        ),
    ]
    fig.legend(
        handles=legend_elements, loc="upper center", ncol=1, bbox_to_anchor=(0.5, 1.02)
    )

    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, "space_time_tradeoff_subplots.png"),
        dpi=300,
        bbox_inches="tight",
    )
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
