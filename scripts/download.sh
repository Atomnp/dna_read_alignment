#!/bin/bash

# ==============================================
# Bash script to download human genome (GRCh38)
# and prepare it for viewing in VS Code
# ==============================================

# Set working directory
WORKDIR="./human_genome"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit

echo "Working directory: $WORKDIR"

# Step 1: Download full genome FASTA (GRCh38)
echo "Downloading human genome (GRCh38) from NCBI..."
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz -O human_genome.fasta.gz

# Step 2: Unzip the FASTA file
echo "Unzipping FASTA file..."
gunzip -f human_genome.fasta.gz

# Step 3: Optional - Rename for convenience
mv -f GCA_000001405.15_GRCh38_genomic.fna human_genome.fasta 2>/dev/null || true

# Step 4: Check the file size
echo "FASTA file size:"
ls -lh human_genome.fasta

# Step 5: Count total bases (optional)
echo "Counting total bases (A/T/C/G)..."
grep -v ">" human_genome.fasta | wc -c

# Step 6: Make a VS Code friendly version (wrap lines at 100 bases)
echo "Wrapping sequence lines at 100 characters for VS Code..."
fold -w 100 human_genome.fasta > human_genome_wrapped.fasta

echo "Done!"
echo "You can now open: $WORKDIR/human_genome_wrapped.fasta in VS Code"
