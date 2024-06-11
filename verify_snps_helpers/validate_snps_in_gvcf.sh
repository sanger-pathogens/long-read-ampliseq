#!/bin/bash

query_file="expected_snps_from_illumina_sequencing"

directory="../results/variants"

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' not found."
    exit 1
fi

# Find and sort the .gz files by their numeric value
files=$(find "$directory" -type f -name "*.gz" | sort -t_ -k2,2n)

# Loop over each sorted .gz file
for file in $files; do
    #echo "Processing $file"

    example_file="${file%.gz}"

    # Unzip the file
    gunzip "$file"

    total=0
    found=0
    missing=0
    pass_count=0

    while IFS= read -r line; do
        ((total++))
        result=$(grep "$line" "$example_file" | grep "PASS")
        if [ -z "$result" ]; then
            #echo "$line: missing"
            ((missing++))
        else
            #echo "$result"
            ((found++))
        fi
    done < $query_file

    # Count total PASS entries in the file
    pass_count=$(grep -c "PASS" "$example_file")

    # Re-gzip the file
    gzip "$example_file"

    #echo "$example_file  total $pass_count"
    echo "$example_file  $found/$total found $pass_count"
done