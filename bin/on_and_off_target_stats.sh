#!/bin/bash

sample_id=$1
on_target_bam=$2
off_target_bam=$3
output_file=$4

echo "Generating on and off-target statistics for $sample_id..."

echo "Counting on-target reads..."
on_target_count=$(samtools view -c $on_target_bam)
echo "Counting off-target reads..."
off_target_count=$(samtools view -c $off_target_bam)
echo "Counting total reads..."
total_count=$((on_target_count + off_target_count))

echo "Calculating on-target percentage..."
on_target_perc=$(awk -v c=$on_target_count -v t=$total_count 'BEGIN {printf "%0.2f", c * 100 / t}')
echo "Calculating off-target percentage..."
off_target_perc=$(awk -v c=$off_target_count -v t=$total_count 'BEGIN {printf "%0.2f", c * 100 / t}')

echo "Outputting to $output_file..."
echo -e "Name,On-target count,Off-target count,On-target percentage,Off-target percentage\\n$sample_id,$on_target_count,$off_target_count,$on_target_perc%,$off_target_perc%" > $output_file
echo "Done!"
