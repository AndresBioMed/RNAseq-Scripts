#!/bin/bash

# Set the directory containing the "*.gz" files
input_folder="/home/andresunix/rnaseq/rawdata/"  # Replace with your actual input folder

# Set the output folder
output_folder="/home/andresunix/rnaseq/kallisto/"

# Set the index file
index_file="/home/andresunix/rnaseq/index/Homo_sapiens.GRCh38.cdna.all.index"

echo "--> multiKallisto will create subdirectories for all your samples while they are being aligned"
echo "Input folder: $input_folder"
echo "Index file: $index_file"
# Set the number of threads
threads=8

# Iterate through all ".gz" files in the input folder
for input_file in "$input_folder"/*.gz; do
    # Extract the base name of the input file (without the path and extension)
    base_name=$(basename -s .fastq.gz "$input_file")

    # Create an output folder for the sample
    sample_output_folder="$output_folder/$base_name"
    mkdir -p "$sample_output_folder"
    echo "-> The sample $base_name is being aligned now by Kallisto"

    # Run kallisto quant for each input file
    kallisto quant -i "$index_file" -o "$sample_output_folder" -t "$threads" --single -l 250 -s 30 "$input_file" > "$sample_output_folder/$base_name.log" 2>&1
    echo "-> Kallisto has finished aligning $base_name, a log file has been also produced"
done

echo "--> multiKallisto has now Finished processing all your samples"