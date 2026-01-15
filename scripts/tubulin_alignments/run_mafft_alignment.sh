#!/bin/bash
# Run MAFFT alignment on beta-tubulin FASTA files
#
# Setup (one-time):
#   conda env create -f envs/mafft_env.yml
#
# Usage:
#   conda activate mafft_env
#   bash scripts/tubulin_alignments/run_mafft_alignment.sh

set -e

INPUT_DIR="data/tubulin_alignments/tubulin_protein_seqs"
OUTPUT_DIR="data/tubulin_alignments/aligned"

# Check mafft is available
if ! command -v mafft &> /dev/null; then
    echo "Error: mafft not found. Please activate the conda environment:"
    echo "  conda activate mafft_env"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each FASTA file
for fasta in "$INPUT_DIR"/*.fa; do
    gene=$(basename "$fasta" .fa)
    output="$OUTPUT_DIR/${gene}_aligned.fa"

    echo "Aligning: $gene"
    mafft --auto "$fasta" > "$output"
done

echo "Alignment complete!"
echo "Output files in: $OUTPUT_DIR"
