RAW_DATADIR = "inputs/raw_data"
SAMPLES = ["2019_CRN1_1_1", "2019_CRN1_1_2", "2019_CRN1_1_3", "2019_CRN1_1_4",
           "2019_CRN1_1_5", "2019_CRN1_2_1", "2019_CRN1_2_2", "2019_CRN1_2_3",
           "2019_CRN1_2_4", "2019_CRN1_2_5", "2019_SNC2_5_1", "2019_SNC2_5_2",
           "2019_SNC2_5_3", "2019_SNC2_5_4", "2019_SNC2_5_5", "2019_SNC2_6_1",
           "2019_SNC2_6_2", "2019_SNC2_6_3", "2019_SNC2_6_4", "2019_SNC2_6_5",
           "2019_RRV1_9_1", "2019_RRV1_9_2", "2019_RRV1_9_3", "2019_RRV1_9_4",
           "2019_RRV1_9_5", "2019_RRV1_10_1", "2019_RRV1_10_2", "2019_RRV1_10_3",
           "2019_RRV1_10_4", "2019_RRV1_10_5", "2019_SNC1_13_1", "2019_SNC1_13_2",
           "2019_SNC1_13_3", "2019_SNC1_13_4", "2019_SNC1_13_5", "2019_SNC1_14_1",
           "2019_SNC1_14_2", "2019_SNC1_14_3", "2019_SNC1_14_4", "2019_SNC1_14_5",
           "2019_AS2_19_1", "2019_AS2_19_2", "2019_AS2_19_3", "2019_AS2_19_4", 
           "2019_AS2_19_5", "2019_AS2_20_1", "2019_AS2_20_2", "2019_AS2_20_3",
           "2019_AS2_20_4", "2019_AS2_20_5", "2019_RRV2_21_1", "2019_RRV2_21_2",
           "2019_RRV2_21_3", "2019_RRV2_21_4", "2019_RRV2_21_5", "2019_RRV2_22_1",
           "2019_RRV2_22_2", "2019_RRV2_22_3", "2019_RRV2_22_4", "2019_RRV2_22_5",
           "2019_AS1_26_1", "2019_AS1_26_2", "2019_AS1_26_3", "2019_AS1_26_4",
           "2019_AS1_26_5", "2019_AS1_27_1", "2019_AS1_27_2", "2019_AS1_27_3",
           "2019_AS1_27_4", "2019_AS1_27_5", "2019_RRV3_29_1", "2019_RRV3_29_2",
           "2019_RRV3_29_3", "2019_RRV3_29_4", "2019_RRV3_29_5", "2019_RRV3_30_1",
           "2019_RRV3_30_2", "2019_RRV3_30_3", "2019_RRV3_30_4", "2019_RRV3_30_5",
           "2019_AV1_33_1", "2019_AV1_33_2", "2019_AV1_33_3", "2019_AV1_33_4",
           "2019_AV1_33_5", "2019_AV1_34_1", "2019_AV1_34_2", "2019_AV1_34_3",
           "2019_AV1_34_4", "2019_AV1_34_5", "2019_AV2_38_1", "2019_AV2_38_2", 
           "2019_AV2_38_3", "2019_AV2_38_4", "2019_AV2_38_5", "2019_AV2_39_1",
           "2019_AV2_39_2", "2019_AV2_39_3", "2019_AV2_39_4", "2019_AV2_39_5",
           "2019_SRH1_41_1", "2019_SRH1_41_2", "2019_SRH1_41_3", "2019_SRH1_41_4",
           "2019_SRH1_41_5", "2019_SRH1_42_1", "2019_SRH1_42_2", "2019_SRH1_42_3",
           "2019_SRH1_42_4", "2019_SRH1_42_5", "2019_SMV2_45_1", "2019_SMV2_45_2",
           "2019_SMV2_45_3", "2019_SMV2_45_4", "2019_SMV2_45_5", "2019_SMV2_48_1",
           "2019_SMV2_48_2", "2019_SMV2_48_3", "2019_SMV2_48_4", "2019_SMV2_48_5",
           "2019_SMV1_50_1", "2019_SMV1_50_2", "2019_SMV1_50_3", "2019_SMV1_50_4",
           "2019_SMV1_50_5", "2019_SMV1_51_1", "2019_SMV1_51_2", "2019_SMV1_51_3",
           "2019_SMV1_51_4", "2019_SMV1_51_5", "2019_OR1_53_1", "2019_OR1_53_2",
           "2019_OR1_53_3", "2019_OR1_53_4", "2019_OR1_53_5", "2019_OR1_54_1",
           "2019_OR1_54_2", "2019_OR1_54_3", "2019_OR1_54_4", "2019_OR1_54_5", 
           "2019_OR2_57_1", "2019_OR2_57_2", "2019_OR2_57_3", "2019_OR2_57_4",
           "2019_OR2_57_5", "2019_OR2_58_1", "2019_OR2_58_2", "2019_OR2_58_3",
           "2019_OR2_58_4", "2019_OR2_58_5", "2019_inoculum_SRH1", "2019_inoculum_AS1",
           "2019_inoculum_SMV2", "2019_inoculum_AS2", "2019_inoculum_RRV1", "2019_inoculum_SMV1"]

rule all:
    input:
       expand("outputs/gather/{sample}_gather.csv", sample = SAMPLES),
       "outputs/counts_all_organisms/raw_counts.tsv",

rule cat_fastq:
    output: 'inputs/cat/{sample}.fq.gz'
    params: indir = RAW_DATADIR
    shell:'''
    cat {params.indir}/{wildcards.sample}_UMI_S*_L00*_R1_001.fastq.gz > {output} 
    '''

rule umi_extract:
    input: "inputs/cat/{sample}.fq.gz"
    output: "outputs/umi_extract/{sample}.fq.gz"
    conda: "envs/umitools.yml"
    shell:'''
    umi_tools extract --stdin={input} --bc-pattern=NNNNNN --log={output}.log --stdout {output}
    '''

rule fastqc:
    input: "outputs/umi_extract/{sample}.fq.gz"
    output: "outputs/fastqc_umi/{sample}.html"
    params: outdir = "outputs/fastqc_umi"
    conda: "envs/fastqc.yml"
    shell:'''
    fastqc -o {params.outdir} -t 8 --nogroup {input}
    '''

rule bbduk_trim:
    input: 
        reads = "outputs/umi_extract/{sample}.fq.gz",
        polya = "inputs/polya.fa",
        adapters = "inputs/truseq_rna.fa.gz"
    output: "outputs/bbduk/{sample}.fq"
    conda: "envs/bbmap.yml"
    shell:'''
    bbduk.sh in={input.reads} out={output} ref={input.polya},{input.adapters} \
    k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
    '''

rule fastqc_trim:
    input: "outputs/bbduk/{sample}.fq"
    output: "outputs/fastqc_trim/{sample}.html"
    params: outdir = "outputs/fastqc_trim"
    conda: "envs/fastqc.yml"
    shell:'''
    fastqc -o {params.outdir} -t 8 --nogroup {input}
    '''

################################################
## Sourmash characterization
################################################

rule sourmash_compute:
    input: "outputs/bbduk/{sample}.fq"
    output: "outputs/sigs/{sample}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compute -k 21,31,51 --track-abundance --scaled 2000 -o {output} {input}
    '''

rule sourmash_compute_scaled100:
    input: "outputs/bbduk/{sample}.fq"
    output: "outputs/sigs_100/{sample}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compute -k 21,31,51 --track-abundance --scaled 100 -o {output} {input}
    '''

rule download_gather_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.tar.gz"
    shell:'''
    wget -O {output} https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k31.tar.gz
    '''

rule untar_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.sbt.json"
    input:  "inputs/gather_databases/genbank-d2-k31.tar.gz"
    params: outdir = "inputs/gather_databases"
    shell: '''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_rna:
    output: "inputs/gather_databases/euk_rna_k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/qk5th/download
    '''

rule untar_rna:
    output: "inputs/gather_databases/euk_rna_k31.sbt.json"
    input:  "inputs/gather_databases/euk_rna_k31.tar.gz"
    params: outdir = "inputs/gather_databases"
    shell: '''
    tar xf {input} -C {params.outdir}
    '''

rule sourmash_gather:
    input: 
        sig="outputs/sigs/{sample}.sig",
        rc212="inputs/gather_databases/rc212.sig",
        db1="inputs/gather_databases/euk_rna_k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
    output:
        csv="outputs/gather/{sample}_gather.csv",
        matches="outputs/gather/{sample}_matches.sig",
        un="outputs/gather/{sample}_un.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.rc212} {input.db1} {input.db2}
    '''


####################################################
## All organism gene counts
####################################################


rule star_index_genome_all:
    input:
        genome = 'inputs/all_organisms/all.fna',
        gtf = 'inputs/all_organisms/all.gtf'
    params: input_dir = 'inputs/all_organisms' 
    conda: 'envs/star.yml'
    output: 'inputs/all_organisms/SAindex'
    shell:'''
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {params.input_dir} \
         --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang  99
    '''

rule star_align_all:
    #v2.5.2a
    input:
        reads = 'outputs/bbduk/{sample}.fq',
        genome_index = 'inputs/all_organisms/SAindex'
    output: 'outputs/star_all_organisms/{sample}Aligned.sortedByCoord.out.bam' 
    params: 
        out_prefix = lambda wildcards: 'outputs/star_all_organisms/' + wildcards.sample,
        genome_dir = 'inputs/all_organisms'
    conda: 'envs/star.yml'
    shell:'''
    STAR --runThreadN 4 --genomeDir {params.genome_dir}      \
        --readFilesIn {input.reads} --outFilterType BySJout  \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8    \
        --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 \
        --alignIntronMax 1000000 --alignMatesGapMax 1000000  \
        --outSAMattributes NH HI NM MD --outSAMtype BAM      \
        SortedByCoordinate --outFileNamePrefix {params.out_prefix}
    '''

rule index_bam_all:
    input: 'outputs/star_all_organisms/{sample}Aligned.sortedByCoord.out.bam'
    output: 'outputs/star_all_organisms/{sample}Aligned.sortedByCoord.out.bam.bai'
    conda: 'envs/samtools.yml'
    shell:'''
    samtools index {input}
    '''
    
rule dedup_UMI_all:
    input: 
        bam = 'outputs/star_all_organisms/{sample}Aligned.sortedByCoord.out.bam',
        bai = 'outputs/star_all_organisms/{sample}Aligned.sortedByCoord.out.bam.bai'
    output: 'outputs/dedup_all_organisms/{sample}.dedup.bam'
    conda: 'envs/umitools.yml'
    shell:'''
    umi_tools dedup -I {input.bam} --output-stats=deduplicated -S {output}
    '''

rule htseq_count_all:
    input:
        bam = 'outputs/dedup_all_organisms/{sample}.dedup.bam',
        gtf = 'inputs/all_organisms/all.gtf'
    output: "outputs/htseq_all_organisms/{sample}_readcounts.txt"
    conda: "envs/htseq.yml"
    shell:'''
    htseq-count -m intersection-nonempty --type gene -s yes -f bam -r pos {input.bam} {input.gtf} > {output}
    '''

rule make_counts_all:
    input: expand("outputs/htseq_all_organisms/{sample}_readcounts.txt", sample = SAMPLES)
    output: "outputs/counts_all_organisms/raw_counts.tsv"
    conda: "envs/tidyverse.yml"
    script: "scripts/make_raw_counts.R"
