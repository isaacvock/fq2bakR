#!/bin/bash

cpus=$1 # number of cpus
sample=$2 # sample names
input=$3 # input fastqs
output=$4 # output
annotation=$5 # gtf file location
format=$6 # Paired or single end

# Define path to sample data and start in correct dir
    cd $MASTER_DIR
    cd "$sample".dir

if ["$format" = "PE"]; then
    # Copy .fastq or .fastq.gz files into dir and create individual r1 and r2 fastq files:
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

            cat *R1* > "$sample"_1.fastq &
            cat *R2* > "$sample"_2.fastq
            wait

            rm -f *R1*
            rm -f *R2*
        else
            echo "* Skipping copy fastq step"
        fi


    # remove duplicate fastq reads:
        if [ "$step_fastuniq" = "TRUE" ]; then
            if [[ -n $FASTUNIQ_MOD ]]; then module load ${FASTUNIQ_MOD}; fi

            $FASTUNIQ \
                -i <(echo "$sample"_1.fastq; echo "$sample"_2.fastq) \
                -o "$sample"_1u.fastq \
                -p "$sample"_2u.fastq &&

             echo "* fastquniq finished for sample " $sample

            in1="$sample"_1u.fastq
            in2="$sample"_2u.fastq
            # Make files to archive:
            rm "$sample"_{1,2}.fastq

        else
            echo "* Skipping fastuniq step"

            in1="$sample"_1.fastq
            in2="$sample"_2.fastq
        fi


    # trim reads
        if [ "$step_cutadapt" = "TRUE" ]; then
            echo "* Running cutadapt in parallel mode for sample $sample"
    

                cutadapt \
                    -a AGATCGGAAGAGC \
                    -A AGATCGGAAGAGC \
                    --minimum-length=20 \
                    --cores="$cpus" \
                    -o "$sample".t.r1.fastq \
                    -p "$sample".t.r2.fastq \
                    "$in1" "$in2" &&
                echo "* cutadapt finished for " $sample



            rm $in1
            rm $in2
        else
            echo "* Skipping cutadapt step"
        fi
elif ["$format" = "SE"]; then

else

fi
