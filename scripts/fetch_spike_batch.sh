#!/bin/bash

mkdir -p data/raw

for id in $(cat data/ids.txt); do
  echo "Fetching $id..."
  if ! bio fetch --id "$id" --db protein --format fasta > "data/raw/${id}.fasta"; then
    echo "Failed to fetch $id"
  else
    echo "Saved: data/raw/${id}.fasta"
  fi
done
#!/bin/bash

mkdir -p data/raw

while read id; do
    echo "Fetching $id..."
    bio fetch --id "$id" --db protein --format fasta >  "data/raw/${id}.fasta"
done < data/ids.txt
