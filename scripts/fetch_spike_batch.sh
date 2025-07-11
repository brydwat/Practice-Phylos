#!/bin/bash

# Path to ids.txt (same folder as this script)
ID_FILE="../data/ids.txt"

# Path to output directory
OUT_DIR="../data/raw"

# Ensure output directory exists (optional if it already exists)
mkdir -p "$OUT_DIR"

# Loop through IDs and fetch
for id in $(cat "$ID_FILE"); do
  echo "Fetching $id..."
  if ! bio fetch --id "$id" --db protein --format fasta > "$OUT_DIR/${id}.fasta"; then
    echo "Failed to fetch $id"
  else
    echo "Saved: $OUT_DIR/${id}.fasta"
  fi
done
