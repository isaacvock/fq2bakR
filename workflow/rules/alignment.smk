if FORMAT == 'PE':


    rule preprocess:
        input:
            get_input_fastqs
        output:
            "results/fastq_cut/{sample}.t.r1.fastq",
            "results/fastq_cut/{sample}.t.r2.fastq"
        log:
            "logs/preprocess/{sample}.log"
        params:
            shellscript=workflow.source_path("../scripts/preprocess_all.sh"),
            format = config["FORMAT"],
            adapter1 = config["adapter1"],
            adapter2 = config["adapter2"]
        threads: workflow.cores
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {input} {params.format} {output} {params.adapter1} {params.adapter2}
            """

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
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h3n} {params.h3n_path} {params.muts} {params.yale} {input} {output}
                """

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
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output}
                """

else:

    rule preprocess:
        input:
            get_input_fastqs
        output:
            "results/fastq_cut/{sample}.t.fastq"
        log:
            "logs/sort_filter/{sample}.log"
        params:
            shellscript=workflow.source_path("../scripts/preprocess_all.sh"),
            format = config["FORMAT"],
            adapter1 = config["adapter1"]
        threads: workflow.cores
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {input} {params.format} {output} {params.adapter1}
            """

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
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h3n} {params.h3n_path} {params.muts} {params.yale} {input} {output}
                """
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
                {params.shellscript} {threads} {wildcards.sample} {params.format} {params.strand} {params.chr} {params.h2} {input} {output}
                """