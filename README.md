# Genome Pattern Matching & Assembly Algorithms

This project implements and benchmarks different string matching algorithms for genome sequence analysis. It includes implementations of Brute Force, Naive Suffix Array, Suffix Array Induced Sorting (SA-IS), and FM-Index algorithms.

## Project Structure

- `src/brute_force.py`: Naive pattern matching implementation.
- `src/suffix_array.py`: Naive Suffix Array construction ($O(N^2 \log N)$) and search.
- `src/sa_is.py`: Linear time Suffix Array construction using SA-IS algorithm ($O(N)$).
- `src/fm_index.py`: FM-Index implementation using SA-IS.
- `src/generate_reads.py`: Script to simulate NGS reads from reference genomes (Batch mode).
- `src/benchmark_algorithms.py`: Script to benchmark time and memory usage of the algorithms (Batch mode).
- `src/plot_benchmark.py`: Script to visualize benchmark results.
- `src/utils.py`: Common utility functions (e.g., FASTA reader).

## Usage

### 1. Generate Simulated Reads

Generate simulated reads for all genomes in the `genome_data` directory. This script creates two sets of reads for each genome:
- **500 reads**: Used for benchmarking the Brute Force algorithm (which is slow).
- **200,000 reads**: Used for benchmarking optimized algorithms (SA-IS, Naive SA, FM-Index).

```bash
python src/generate_reads.py --input_dir ./genome_data --output_dir ./simulated_reads
```

**Arguments:**
- `--input_dir`: Directory containing reference genome FASTA files (default: `./genome_data`).
- `--output_dir`: Directory to save generated read files (default: `./simulated_reads`).

### 2. Run Benchmarks

Run the benchmark script to compare the performance (Time & Memory) of the algorithms across all genomes. This will execute the algorithms against the generated reads and save the results to a CSV file.

```bash
python src/benchmark_algorithms.py --genome_dir ./genome_data --reads_dir ./simulated_reads --output_csv ./benchmark_results.csv
```

**Arguments:**
- `--genome_dir`: Directory containing reference genomes.
- `--reads_dir`: Directory containing the simulated reads generated in step 1.
- `--output_csv`: Path to save the benchmark results CSV.
- `--measure-memory`: (Optional) Enable memory usage tracking (slower).

### 3. Generate Plots

Visualize the benchmark results using the plotting script. This will generate academic-style charts showing indexing time, search time complexity, and total execution time.

```bash
python src/plot_benchmark.py --csv ./benchmark_results.csv --output ./plots
```

**Arguments:**
- `--csv`: Path to the benchmark results CSV file.
- `--output`: Directory to save the generated plots.

## Algorithms Implemented

1.  **Brute Force**: Naive string matching.
2.  **Naive Suffix Array**: $O(N^2 \log N)$ construction.
3.  **SA-IS**: Linear time $O(N)$ Suffix Array construction.
4.  **FM-Index**: Compressed full-text index based on the Burrows-Wheeler Transform.
