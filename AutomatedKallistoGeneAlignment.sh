#!/bin/bash

conda activate rnaseq

echo "----Wellcome to Automated Kallisto Gene alignment [AKG]----"
echo "This script will align all your samples in a new folder using Kallisto, and checking their quality with FastQC and MultiQC."

mkdir -p /home/andresunix/rnaseq/new_AKG
mkdir -p /home/andresunix/rnaseq/new_AKG/fastqc
mkdir -p /home/andresunix/rnaseq/new_AKG/kallisto
mkdir -p /home/andresunix/rnaseq/new_AKG/index

# Prompt the user for the input folder
read -p "Enter the absolute path to the input folder (containing *.gz files): " input_folder

# Prompt the user for the output folder
read -p "Enter the absolute path to the reference genome (containing a *.fa file): " genome_file

echo "Input folder: $input_folder"
echo "Genome file: $genome_file"

# Analyze .gz files with FastQC
echo "--> AKG will now analyze the quality of your samples with FastQC"
cd $input_folder
fastqc *.gz -t 8
mv *fastqc* /home/andresunix/rnaseq/new_AKG/fastqc
echo "--> AKG has finished analyzing the quality of your samples with FastQC"
echo "--> AKG will now create an index based on your refernce genome"
# eate the index for the reference genome
cd /home/andresunix/rnaseq/new_AKG/index
kallisto index -i Homo_sapiens.GRCh38.cdna.all.index $genome_file/Homo_sapiens.GRCh38.cdna.all.fa 
echo "--> AKG has finished the index"
# Set the number of threads
threads=8

# Iterate through all ".gz" files in the input folder
for input_file in "$input_folder"/*.gz; do
    # Extract the base name of the input file (without the path and extension)
    base_name=$(basename -s .fastq.gz "$input_file")

    # Create an output folder for the sample
    sample_output_folder="/home/andresunix/rnaseq/new_AKG/kallisto/$base_name"
    mkdir -p "$sample_output_folder"
    echo "-> The sample $base_name is being aligned now by Kallisto"

    # Run kallisto quant for each input file
    kallisto quant -i "/home/andresunix/rnaseq/new_AKG/index/Homo_sapiens.GRCh38.cdna.all.index" -o "$sample_output_folder" -t "$threads" --single -l 250 -s 30 "$input_file" > "$sample_output_folder/$base_name.log" 2>&1
    echo "-> Kallisto has finished aligning $base_name, a log file has also been produced"
done

echo "--> AKG has now Finished processing all your samples"
echo "Summarising results via MultiQ"

cd /home/andresunix/rnaseq/new_AKG
multiqc -d .
echo "AKG has finished, the final report has ben produced alongside with the pseudoalignments"
