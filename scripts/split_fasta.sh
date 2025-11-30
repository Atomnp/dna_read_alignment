#!/bin/bash

# ==============================================
# Bash script to split a large FASTA file
# into smaller files for easy viewing
# ==============================================

# Input FASTA file
INPUT_FASTA="./human_genome/human_genome.fasta"

# Output directory
OUTDIR="./human_genome/fasta_chunks"
mkdir -p "$OUTDIR"

# Number of chunks you want
NUM_CHUNKS=20

# Count total lines (excluding header lines starting with '>')
TOTAL_LINES=$(grep -v "^>" "$INPUT_FASTA" | wc -l)
echo "Total sequence lines (excluding headers): $TOTAL_LINES"

# Calculate lines per chunk
LINES_PER_CHUNK=$(( (TOTAL_LINES + NUM_CHUNKS - 1) / NUM_CHUNKS ))
echo "Lines per chunk: $LINES_PER_CHUNK"

# Temporary file to hold only sequences (no headers)
grep -v "^>" "$INPUT_FASTA" > "$OUTDIR/temp_sequence.txt"

# Split the sequence into chunks
split -l $LINES_PER_CHUNK "$OUTDIR/temp_sequence.txt" "$OUTDIR/chunk_"

# Add back the header for each chunk
cd "$OUTDIR" || exit
for f in chunk_*; do
    HEADER=">partial_sequence_${f}"
    mv "$f" "${f}.fasta"
    sed -i "1i$HEADER" "${f}.fasta"
done

# Remove temporary file
rm temp_sequence.txt

echo "Done! Files saved in $OUTDIR"
echo "You can now open these smaller FASTA files in VS Code"
