#!/bin/bash
# Main sript for sample preprocessing

# This allows to be script run as: sbatch sample_muts.sh configFile sampleName
    if [[ -n $1 ]]; then config=$1; fi
    if [[ -n $2 ]]; then sample=$2; fi

# Source the paths and variables:
    source $config

    cd $MASTER_DIR

    if  [ "$step_copyfq" = "TRUE" ] || \
        [ "$step_fastuniq" = "TRUE" ] || \
        [ "$step_cutadapt" = "TRUE" ]; then

        $PREPROCESS_SCRIPT $config $sample &&
        echo "** Preprocessing finished for sample $sample" ||
        ( echo "!!! Error in preprocessing"; exit 1 )
    else
        echo "** Skipping preprocessing for sample $sample"
    fi


    if [ "$step_align" = "TRUE" ]; then
        $ALIGNMENT_SCRIPT $config $sample &&
        echo "** Alignment completed for sample $sample" ||
        ( echo "!!! Error in aligning"; exit 1 )
    else
        echo "** Skipping alignment for sample $sample"
    fi


    if [ "$step_readfilter" = "TRUE" ]; then
        $SORT_FILTER_SCRIPT $config $sample &&
        echo "** Aligned reads filtered for sample $sample" ||
        ( echo "!!! Error in filtering"; exit 1 )
    else
        echo "** Skipping aligned reads filtering for sample $sample"
    fi


    if [ "$step_feature_count" = "TRUE" ]; then
        $FEATURE_COUNT_SCRIPT $config $sample &&
        echo "** HTSeq feature counting finished for sample $sample"||
        ( echo "!!! Error in feature counting"; exit 1 )
     else
        echo "** Skipping HTseq counting for sample $sample"
    fi


    if [ "$step_fragment" = "TRUE" ]; then
        $FRAGMENT_SCRIPT $config $sample &&
        echo "** Aligned .sam file fragmented for sample $sample" ||
        ( echo "!!! Error in fragmentation"; exit 1 )
     else
        echo "** Skipping .sam fragmentation for sample $sample"
    fi
