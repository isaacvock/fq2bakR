if FORMAT == 'PE':

    rule cutadapt:
        input:
            get_input_fastqs
        output:
            fastq1="results/fastq_cut/{sample}.t.r1.fastq",
            fastq2="results/fastq_cut/{sample}.t.r2.fastq",
            qc = "results/fastq_cut/{sample}.qc.txt",
        params:
            adapters=config["adapter"],
            extra=config["extra"], 
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8  # set desired number of threads here
        wrapper:
            "v1.25.0/bio/cutadapt/pe"       

    if ALIGNER:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.r1.fastq",
                "results/fastq_cut/{sample}.t.r2.fastq"
            output:
                "results/bams/{sample}Aligned.out.bam",
                "results/bams/{sample}.sam"
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
                    str(config["annotation"]), config["star_params"]
                ),
            conda:
                "../envs/star.yaml"
            threads: 24
            script: 
                "../scripts/star-align.py"

    else:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.r1.fastq",
                "results/fastq_cut/{sample}.t.r2.fastq"
            output:
                "results/bams/{sample}Aligned.out.bam",
                "results/bams/{sample}.sam"
            log:
                "logs/align/{sample}.log"
            params:
                shellscript = workflow.source_path("../scripts/hisat2.sh"),
                format = config["FORMAT"],
                strand = config["strandedness"],
                chr = config["chr_tag"],
                h2 = config["HISAT2"]
            threads: workflow.cores
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output} 1> {log} 2>&1
                """

else:

    rule cutadapt:
        input:
            get_input_fastqs
        output:
            fastq="results/fastq_cut/{sample}.t.fastq",
            qc = "results/fastq_cut/{sample}.qc.txt",
        params:
            adapters=config["adapter"],
            extra=config["extra"], 
        log:
            "logs/cutadapt/{sample}.log",
        threads: 8  # set desired number of threads here
        wrapper:
            "v1.25.0/bio/cutadapt/se" 

    if ALIGNER:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.fastq",
            output:
                "results/bams/{sample}Aligned.out.bam",
                "results/bams/{sample}.sam"
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
                    str(config["annotation"]), config["star_params"]
                ),
            conda:
                "../envs/star.yaml"
            threads: 24
            script: 
                "../scripts/star-align.py"

    else:
        rule align:
            input:
                "results/fastq_cut/{sample}.t.fastq",
            output:
                "results/bams/{sample}Aligned.out.bam",
                "results/bams/{sample}.sam"
            log:
                "logs/align/{sample}.log"
            params:
                shellscript = workflow.source_path("../scripts/hisat2.sh"),
                format = config["FORMAT"],
                strand = config["strandedness"],
                chr = config["chr_tag"],
                h2 = config["HISAT2"]
            threads: workflow.cores
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output} 1> {log} 2>&1
                """
