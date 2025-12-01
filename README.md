# Genome Pattern Matching & Assembly Algorithms

This project implements and benchmarks different string matching algorithms for genome sequence analysis. It includes implementations of Brute Force, Naive Suffix Array, and Suffix Array Induced Sorting (SA-IS) algorithms.

## Project Structure

- `src/brute_force.py`: Naive pattern matching implementation.
- `src/suffix_array.py`: Naive Suffix Array construction ($O(N^2 \log N)$) and search.
- `src/sa_is.py`: Linear time Suffix Array construction using SA-IS algorithm ($O(N)$).
- `src/generate_reads.py`: Script to simulate NGS reads from a reference genome.
- `src/benchmark_algorithms.py`: Script to benchmark time and memory usage of the algorithms.
- `src/utils.py`: Common utility functions (e.g., FASTA reader).

## Usage

### 1. Generate Simulated Reads

First, generate a set of simulated reads from a reference genome. This mimics the output of a DNA sequencer.

```bash
python src/generate_reads.py --fasta <path_to_genome> --num_reads <number> --read_len <length> --output <output_path>
```

**Arguments:**
- `--fasta`: Path to the input reference genome FASTA file.
- `--num_reads`: Number of reads to generate.
- `--read_len`: Length of each read (default: 100).
- `--output`: Path where the generated reads FASTA file will be saved.

### 2. Run Benchmarks

Run the benchmark script to compare the performance (Time & Memory) of the implemented algorithms in mapping the generated reads back to the reference genome.

```bash
python src/benchmark_algorithms.py --fasta <path_to_genome> --reads <path_to_reads> [--skip-brute]
```

**Arguments:**
- `--fasta`: Path to the reference genome FASTA file.
- `--reads`: Path to the reads FASTA file generated in step 1.
- `--skip-brute`: (Optional) Skip the Brute Force algorithm. Recommended for large genomes as it is significantly slower.

## Example: E. coli K-12

Here is a complete example using the *E. coli* K-12 genome.

**Step 1: Generate 1,000 reads of length 50**

```bash
python src/generate_reads.py --fasta ./genome_data/e_coli_k12.fasta --num_reads 1000 --read_len 50 --output ./genome_data/e_coli_reads.fasta
```

**Step 2: Run the benchmark**

```bash
python src/benchmark_algorithms.py --fasta ./genome_data/e_coli_k12.fasta --reads ./genome_data/e_coli_reads.fasta
```

*Note: If the benchmark takes too long, you can skip the brute force method:*

```bash
python src/benchmark_algorithms.py --fasta ./genome_data/e_coli_k12.fasta --reads ./genome_data/e_coli_reads.fasta --skip-brute
```
