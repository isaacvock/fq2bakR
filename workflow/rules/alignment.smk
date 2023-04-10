## Create STAR index if not available
if config['build_star']:
    rule star_index:
        input:
            fasta=config["genome_fasta"],
            annotation=config["annotation"],
        output:
            directory(config['STAR_index']),
        threads: 4
        params:
            extra="--sjdbGTFfile {} --sjdbOverhang 100".format(str(config["annotation"])),
        log:
            "logs/star_index_genome.log",
        wrapper:
            "v1.25.0/bio/star/index"


## Paired-end reads
if FORMAT == 'PE':

    # Run cutadapt
    rule cutadapt:
        input:
            get_input_fastqs
        output:
            fastq1="results/fastq_cut/{sample}.t.r1.fastq",
            fastq2="results/fastq_cut/{sample}.t.r2.fastq",
            qc = "results/fastq_cut/{sample}.qc.txt",
        params:
            adapters=config["adapter"],
            extra=config["cutadapt_extra"], 
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8  # set desired number of threads here
        wrapper:
            "v1.25.0/bio/cutadapt/pe"       

    # Run hisat-3n
    if ALIGNER:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.r1.fastq",
                "results/fastq_cut/{sample}.t.r2.fastq"
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/align/{sample}.log"
            params:
                shellscript = workflow.source_path("../scripts/hisat_3n.sh"),
                format = config["FORMAT"],
                strand = config["strandedness"],
                chr = config["chr_tag"],
                h3n = config["HISAT_3N"],
                h3n_path = config["hisat3n_path"],
                muts = config["mut_tracks"],
                yale = config["Yale"]
            threads: workflow.cores
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h3n} {params.h3n_path} {params.muts} {params.yale} {input} {output} 1> {log} 2>&1
                """

    # Run STAR
    elif STAR:

        rule align:
            input:
                fq1 = "results/fastq_cut/{sample}.t.r1.fastq",
                fq2 = "results/fastq_cut/{sample}.t.r2.fastq",
                index = config['STAR_index'],
            output:
                aln="results/bams/{sample}Aligned.out.bam",
                reads_per_gene="results/bams/{sample}-ReadsPerGene.out.tab",
                aln_tx="results/bams/{sample}-Aligned.toTranscriptome.out.bam",
            log:
                "logs/bams/{sample}.log",
            params:
                idx=lambda wc, input: input.index,
                extra="--outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 23 --outSAMattributes NH HI AS NM MD --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile {} {}".format(
                    str(config["annotation"]), config["star_extra"]
                ),
            conda:
                "../envs/star.yaml"
            threads: 24
            script: 
                "../scripts/star-align.py"

    # Run hisat2
    else:
        ## Old way of running hisat2 with custom script
       # rule align:
           # input:
           #     "results/fastq_cut/{sample}.t.r1.fastq",
           #     "results/fastq_cut/{sample}.t.r2.fastq"
           # output:
           #     "results/bams/{sample}Aligned.out.bam",
           # log:
           #     "logs/align/{sample}.log"
           # params:
           #     shellscript = workflow.source_path("../scripts/hisat2.sh"),
           #     format = config["FORMAT"],
           #     strand = config["strandedness"],
           #     chr = config["chr_tag"],
           #     h2 = config["HISAT2"]
           # threads: workflow.cores
           # conda:
           #     "../envs/full.yaml"
           # shell:
           #     """
           #     chmod +x {params.shellscript}
           #     {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output} 1> {log} 2>&1
           #     """
        
        # Running hisat2 with Snakemake wrapper
        rule align:
            input:
                reads=["results/fastq_cut/{sample}.t.r1.fastq","results/fastq_cut/{sample}.t.r2.fastq"],
                idx=config["HISAT2"],
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/hisat2_align/{sample}.log",
            params:
                extra=config["hisat2_extra"],
            threads: workflow.cores
            wrapper:
                "v1.25.0/bio/hisat2/align"

## Single end data
else:

    # Run cutadapt
    rule cutadapt:
        input:
            get_input_fastqs
        output:
            fastq="results/fastq_cut/{sample}.t.fastq",
            qc = "results/fastq_cut/{sample}.qc.txt",
        params:
            adapters=config["adapter"],
            extra=config["cutadapt_extra"], 
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8  # set desired number of threads here
        wrapper:
            "v1.25.0/bio/cutadapt/se" 

    # Run hisat-3n
    if ALIGNER:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.fastq",
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/align/{sample}.log"
            params:
                shellscript = workflow.source_path("../scripts/hisat_3n.sh"),
                format = config["FORMAT"],
                strand = config["strandedness"],
                chr = config["chr_tag"],
                h3n = config["HISAT_3N"],
                h3n_path = config["hisat3n_path"],
                muts = config["mut_tracks"],
                yale = config["Yale"]
            threads: workflow.cores
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h3n} {params.h3n_path} {params.muts} {params.yale} {input} {output} 1> {log} 2>&1
                """

    # Run STAR
    elif STAR:
        rule align:
            input:
                fq1 = "results/fastq_cut/{sample}.t.fastq",
                index = config['STAR_index'],
            output:
                aln="results/bams/{sample}Aligned.out.bam",
                reads_per_gene="results/bams/{sample}-ReadsPerGene.out.tab",
                aln_tx="results/bams/{sample}-Aligned.toTranscriptome.out.bam",
            log:
                "logs/bams/{sample}.log",
            params:
                idx=lambda wc, input: input.index,
                extra="--outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 23 --outSAMattributes NH HI AS NM MD --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile {} {}".format(
                    str(config["annotation"]), config["star_extra"]
                ),
            conda:
                "../envs/star.yaml"
            threads: 36
            script: 
                "../scripts/star-align.py"

    # Run hisat2
    else:
        ## Old way of running hisat2 with custom script
        #rule align:
            #input:
            #    "results/fastq_cut/{sample}.t.fastq",
            #output:
            #    "results/bams/{sample}Aligned.out.bam",
            #log:
            #    "logs/align/{sample}.log"
            #params:
            #    shellscript = workflow.source_path("../scripts/hisat2.sh"),
            #    format = config["FORMAT"],
            #    strand = config["strandedness"],
            #    chr = config["chr_tag"],
            #    h2 = config["HISAT2"]
            #threads: workflow.cores
            #conda:
            #    "../envs/full.yaml"
            #shell:
            #    """
            #    chmod +x {params.shellscript}
            #    {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output} 1> {log} 2>&1
            #    """
        
        # Run hisat2 with Snakemake wrapper
        rule align:
            input:
                reads=["results/fastq_cut/{sample}.t.fastq"],
                idx=config["HISAT2"],
            output:
                "results/bams/{sample}Aligned.out.bam",
            log:
                "logs/hisat2_align/{sample}.log",
            params:
                extra=config["hisat2_extra"],
            threads: workflow.cores
            wrapper:
                "v1.25.0/bio/hisat2/align"
