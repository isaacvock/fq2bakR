#!/bin/bash

## Maybe I can have input be folders containing fastq files (of either SE or PE variety).
## Throw out step where I check if gz or not (force them to be fastq files)
## But other problem I have is that for whatever reasons, files get renamed
    ## I guess this is needed to ensure consistency of naming

## Realized I can get file names from fastq folder with some bash
    ## names=($(ls $fastq_dir)); names[0] is one of the read pair, names[1] is other if it exists
    cpus=$1 # number of cpus
    sample=$2 # sample names
    fastq_dir=$3 # input fastqs
    format=$4 # Paired or single end



if ["$format" = "PE"]; then

    output1=$5 # output1
    output2=$6

    # Copy .fastq or .fastq.gz files into dir and create individual r1 and r2 fastq files:

    if [[ -f $(echo ${fastq_dir}/*.fastq | cut -f 1 -d " ") ]]; then
        suffix="fastq"
    elif [[ -f $(echo ${fastq_dir}/*.gz | cut -f 1 -d " ") ]]; then
        suffix="gz"
    else
        echo "!!! No fastq files found at location ${fastq_dir}"
        exit 1
    fi

    fastqs=($(ls ${fastq_dir}))

    # remove duplicate fastq reads:
            fastuniq \
                -i <(echo "${fastqs[0]}"; echo "${fastqs[1]}") \
                -o "$sample"_1u.fastq \
                -p "$sample"_2u.fastq &&

             echo "* fastquniq finished for sample " $sample



    # trim reads
            echo "* Running cutadapt in parallel mode for sample $sample"


                cutadapt \
                    -a AGATCGGAAGAGC \
                    -A AGATCGGAAGAGC \
                    --minimum-length=20 \
                    --cores="$cpus" \
                    -o "$output1" \
                    -p "$output2"\
                    "$sample"_1u.fastq "$sample"_2u.fastq &&
                echo "* cutadapt finished for " $sample



            rm "$sample"_1u.fastq
            rm "$sample"_2u.fastq


elif ["$format" = "SE"]; then

    output1=$5 # output1

    # Copy .fastq or .fastq.gz files into dir and create individual r1:
        if [ "$step_copyfq" = "TRUE" ]; then
            # Check that files exist
            if [[ -f $(echo ${LINK_BASE}/${prefix}${sample}/*.fastq | cut -f 1 -d " ") ]]; then
                suffix="fastq"
            elif [[ -f $(echo ${LINK_BASE}/${prefix}${sample}/*.gz | cut -f 1 -d " ") ]]; then
                suffix="gz"
            else
                echo "!!! No fastq files found at location ${LINK_BASE}/${prefix}${sample}"
                exit 1
            fi


            # parallel copying of files
            if [[ -n $GNUPARALLEL_MOD ]]; then module load ${GNUPARALLEL_MOD}; fi

            parallel -j $cpus "cp {1} ./" :::  ${LINK_BASE}/${prefix}${sample}/*.${suffix}

            echo "* .${suffix} files copied for sample $sample"


            # Parallel decompression
            if [[ $suffix == "gz" ]]; then
                if [[ -n $PIGZ_MOD ]]; then module load ${PIGZ_MOD}; fi

                $PIGZ -d -p "$cpus" *.gz &&
                echo "* .fastq files decompressed for sample $sample"
            fi

            cat *R1* > "$sample"_1.fastq

            rm -f *R1*

        else
            echo "* Skipping copy fastq step"
        fi


        echo "* Skipping fastuniq step for SE data"
        in="$sample"_1.fastq


        if [ "$step_cutadapt" = "TRUE" ]; then
            echo "* Running cutadapt in parallel mode for sample $sample"

            if [[ -n $CUTADAPT_MOD ]]; then module load ${CUTADAPT_MOD}; fi

            if [ "$STL" = "TRUE" ]; then
                $CUTADAPT \
                    -a ^CGATC...GATCGGAAGAGC \
                    -a ^CGATC...GGAAGAGCACAC \
                    -n 2 \
                    --minimum-length=20 \
                    --cores="$cpus" \
                    -o "$sample".t.fastq \
                    "$in" &&
                echo "* cutadapt finished for " $sample
            else
                $CUTADAPT \
                    -a AGATCGGAAGAGC \
                    --minimum-length=20 \
                    --cores="$cpus" \
                    -o "$sample".t.fastq \
                    "$in" &&
                echo "* cutadapt finished for " $sample
            fi


            rm $in
        else
            echo "* Skipping cutadapt step"
        fi

else
    echo "!!! format must be PE or SE !!!"
fi
