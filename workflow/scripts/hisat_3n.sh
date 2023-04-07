#!/bin/bash

# Source the paths and variables:
cpus=$1 # number of cpus
sample=$2 # sample names
format=$3 # Paired or single end
reads=$4 # FR or RF
chr_tag=$5
HISAT_3N_INDEX=$6
HISAT_3N_PATH=$7
mut_tracks=$8

if [ "$format" = "PE" ]; then
    input1=$9
    input2=${10}
    output=${11}
    output2=${12}

else
    input=$9
    output=${10}
    outpu2=${11}
fi



    echo "* Aligning reads with HISAT-3n for sample $sample"

    if [[ "$chr_tag" = "TRUE" ]]; then
        echo "* chr tag will be included in alignment for " $sample
    fi


    if [[ "$format" = "PE" ]]; then
        if [[ "$reads" = "F"]]; then
        
            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness FR \
                -S "$sample".sam
        else
            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -1 "$input1" \
                -2 "$input2" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness RF \
                -S "$sample".sam
    elif [[ "$reads" = "SE" ]]; then
        if [[ "$reads" = "F"]]; then

            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -U "$input1" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness F \
                -S "$sample".sam
        else

            $HISAT_3N_PATH \
                -p "$cpus" \
                -x $HISAT_3N_INDEX \
                -U "$input1" \
                $( if [ $mut_tracks = GA ]; then echo "--base-change G,A"; else echo "--base-change T,C"; fi ) \
                $( if [ $chr_tag = 'TRUE' ]; then echo "--add-chrname "; fi ) \
                --rna-strandness R \
                -S "$sample".sam
    else
        echo "! No FR/RF/F method recognized for " $sample
        exit 1
    fi &&


    samtools view -@ "$cpus" -o "$output" "$sample".sam

    rm "$sample".sam

    echo "* Alignment script finished for " $sample
