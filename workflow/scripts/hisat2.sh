#!/bin/bash

# Source the paths and variables:

cpus=$1 # number of cpus
sample=$2 # sample names
format=$3 # Paired or single end
reads=$4 # F or R
chr_tag=$5
HISAT_INDEX=$6

if [ "$format" = "PE" ]; then
    input1=$7
    input2=$8
    output=$9
    outpu2=${10}

else
    input=$7
    output=$8
    output2=$9
fi


# align trimmed reads to the  genome
    echo "* Aligning reads with HISAT2 for sample $sample"

    if [[ "$chr_tag" = "TRUE" ]]; then
        echo "* chr tag will be included in alignment for " $sample
    fi


    if [[ "$format" = "PE" ]]; then
        if [[ "$reads" = "F"]]; then
        
            hisat2 \
                -p "$cpus" \
                -x $HISAT_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness FR \
                --mp 4,2 \
                -S "$output2"
        else
            hisat2 \
                -p "$cpus" \
                -x $HISAT_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness FR \
                --mp 4,2 \
                -S "$output2"

    elif [[ "$format" = "SE" ]]; then
        if [[ "$reads" = "F"]]; then

            hisat2 \
                -p "$cpus" \
                -x $HISAT_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness FR \
                --mp 4,2 \
                -S "$output2"
        else

            hisat2 \
                -p "$cpus" \
                -x $HISAT_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness FR \
                --mp 4,2 \
                -S "$output2"
    else
        echo "! No PE/SE FORMAT parameter recognized for " $sample
        exit 1
    fi &&

    samtools view -@ "$cpus" -o "$output" "$output2"

    echo "* Alignment script finished for " $sample
