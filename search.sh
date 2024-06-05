#!/bin/bash

query_file="snps"

directory="./results/variants"

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' not found."
    exit 1
fi

# Loop over each .gz file
for file in $(find "$directory" -type f -name "*.gz"); do
    echo "Processing $file"

    example_file="${file%.gz}"

    # Unzip the query file
    gunzip "$example_file.gz"

    total=0
    found=0
    missing=0

    while IFS= read -r line
    do
        ((total++))
        result=$(grep "$line" "$example_file" | grep "PASS")
        if [ -z "$result" ]; then
            #echo "$line: missing"
            ((missing++))
        else
            echo "$result"
            ((found++))
        fi
    done < $query_file

    # Re-gzip the query file
    gzip "$example_file"

    echo "$found/$total found"
done