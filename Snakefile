from scripts.split_fasta_regions import split_fasta
from snakemake.utils import R
from collections.abc import Iterable

import getpass
import os
import glob

# Corrected localrules assignment
localrules: all

REFERENCE = "./reference/GCF_002443255.1_Vdes_3.0_genomic.fna"
RAW_READS_DIR = "/data/lyc_vcf_new/input"
OUTPUT_DIR = "/data/lyc_vcf_next/output_highqc1120"
SCRATCH = "/data/lyc_vcf_next/temp_high"

VDESRef = REFERENCE
SPLITS = range(300)
REGIONS = split_fasta(VDESRef, len(SPLITS))

# Corrected region processing
for region in REGIONS:
    for idx, i in enumerate(REGIONS[region]):
        REGIONS[region][idx] = " -r " + str(i)

SAMPLES = []
for fq1 in glob.glob(f"{RAW_READS_DIR}/*_1.fastq.gz"):
    base_name = os.path.basename(fq1).rsplit("_", 1)[0]
    fq2 = fq1.replace("_1.fastq.gz", "_2.fastq.gz")
    if os.path.exists(fq2):
        SAMPLES.append(base_name)
    else:
        print(f"Warning：sample {base_name} lack the paired _2 file，skip")

SAMPLES = sorted(list(set(SAMPLES)))
print(f"sample_list：{SAMPLES}")


rule all:
    input:
        expand(f"/data/lyc_vcf_next/output/alignments/ngm/{{sample}}.bam", sample=SAMPLES),
        expand(f"/data/lyc_vcf_next/output/alignments/ngm/{{sample}}.bam.bai", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/var/ngm/fix_freebayes_filtered.vcf")
        #expand(f"{OUTPUT_DIR}/meta/{{sample}}.txt",sample=SAMPLES)


rule align_and_dedup:
    input:
        read1="/data/lyc_vcf_new/input/{sample}_1.fastq.gz",
        read2="/data/lyc_vcf_new/input/{sample}_2.fastq.gz",
        ref="/home/lyh/lyc/project_vcf_M4_nextgenmap/reference/GCF_002443255.1_Vdes_3.0_genomic.fna"
    output:
        bam=temp("/data/lyc_vcf_next/output/alignments/ngm/{sample}.dedup.bam"),
        bai=temp("/data/lyc_vcf_next/output/alignments/ngm/{sample}.dedup.bam.bai")
    shell:
        """
        mkdir -p /data/lyc_vcf_next/output_highqc/alignments/ngm/
        mkdir -p /data/lyc_vcf_next/temp/ngm/
        
        echo "___ngm___" 
        ngm -t 8 --qry1 {input.read1} --qry2 {input.read2} --paired -r {input.ref} --local --very-sensitive --rg-id {wildcards.sample} --rg-sm {wildcards.sample} --rg-pl ILLUMINA --rg-lb NEXTERA --rg-cn ZJU | 
        samtools view -Su - | 
        samtools fixmate -u -m - - |\
        samtools sort - -m 40G -T /data/lyc_vcf_next/temp/ngm/{wildcards.sample} -o /data/lyc_vcf_next/temp/ngm/{wildcards.sample}.sorted.bam
        
        samtools markdup -r /data/lyc_vcf_next/temp/ngm/{wildcards.sample}.sorted.bam {output.bam}
        samtools index {output.bam}
        
        samtools quickcheck {output.bam} || (echo "Fail"; exit 1)
        """

rule run_variant:
    input:
        bam=rules.align_and_dedup.output.bam,
        bai=rules.align_and_dedup.output.bai
    output:
        bam="/data/lyc_vcf_next/output/alignments/ngm/{sample}.bam",
        bai="/data/lyc_vcf_next/output/alignments/ngm/{sample}.bam.bai",
        log="/data/lyc_vcf_next/output/alignments/ngm/{sample}.variant.log"
    shell:
        """
        mkdir -p /data/lyc_vcf_next/output/alignments/ngm/
        
        variant {input.bam} -m 50 --bam -o {output.bam} 2> {output.log}
        
        samtools index {output.bam}
        """



rule statsbam:
    input:
        alignment = f"/data/lyc_vcf_next/output/alignments/ngm/{{sample}}.bam"
    output:
        f"/data/lyc_vcf_next/output/meta/{{sample}}.txt"
    threads: 16

    shell:
        """
        mkdir -p /data/lyc_vcf_next/output/meta
		echo {wildcards.sample} > {output}
		samtools depth -a {input.alignment} | awk '{{sum+=$3}} END {{ print "Mean Average Coverage on all sites = ",sum/NR}}' >> {output}
		samtools depth -a -r NW_019211454.1 -r NW_019211455.1 -r NW_019211456.1 -r NW_019211457.1 -r NW_019211458.1 -r NW_019211459.1 -r NW_019211460.1 {input.alignment} | awk '{{sum+=$3}} END {{ print "Mean Average Coverage on 7 chromosomes = ",sum/NR}}' >> {output}
		samtools depth -a -r NC_004454.2 {input.alignment} | awk '{{sum+=$3}} END {{ print "Mean Average Coverage mtDNA = ",sum/NR}}' >> {output}
		samtools flagstat {input.alignment} >> {output}

		"""

rule freebayes:
    input:
        expand(f"/data/lyc_vcf_next/output/alignments/ngm/{{sample}}.bam", sample=SAMPLES)
    output:
        temp(f"{OUTPUT_DIR}/var/ngm/snponly_split/fix_freebayes.{{region}}.vcf")
    params:
        span = lambda wildcards: REGIONS[wildcards.region],
        bams = lambda wildcards, input: " ".join(input),
        filtering = "--min-alternate-count 2 --min-alternate-fraction 0.2 --min-mapping-quality 8 --min-base-quality 5 --use-best-n-alleles 4"
    conda:
        "envs/freebayes_env.yaml"
    shell:
        """
        freebayes {params.filtering} --fasta-reference {VDESRef} {params.span} {params.bams} > {output}
        """

rule merge_vcf:
    input:
        expand(f"{OUTPUT_DIR}/var/ngm/snponly_split/fix_freebayes.{{region}}.vcf", region=REGIONS.keys())
    output:
        f"{OUTPUT_DIR}/var/ngm/fix_snponly_freebayes.vcf"
    conda:
        "envs/merge_vcf_env.yaml"
    shell:
        """
        mkdir -p {OUTPUT_DIR}/var/ngm
        bcftools concat -Oz -o {output}.tmp {input}
        bcftools view -Ov {output}.tmp > {output}
        rm {output}.tmp
        """

rule filter_vcf:
    input:
        vcf = f"{OUTPUT_DIR}/var/ngm/fix_snponly_freebayes.vcf"  
    output:
        vcf_filtered = f"{OUTPUT_DIR}/var/ngm/fix_freebayes_filtered.vcf" 
    params:
        span = "--chr NW_019211454.1 --chr NW_019211455.1 --chr NW_019211456.1 --chr NW_019211457.1 --chr NW_019211458.1 --chr NW_019211459.1 --chr NW_019211460.1",
        filters = "--max-alleles 2 --minQ 40 --minDP 2 --maxDP 50 --max-missing 0.5 --maf 0.05"
    shell:
        """
        echo "start generating VCF file..."
        OUTPUT_FILE="{output.vcf_filtered}"
        OUTPUT_PREFIX="${{OUTPUT_FILE%.vcf}}"
        vcftools --vcf '{input.vcf}' {params.span} {params.filters} --recode --recode-INFO-all --out $OUTPUT_PREFIX 2> $OUTPUT_PREFIX.log
        echo "finish generating VCF file..."
        mv $OUTPUT_PREFIX.recode.vcf {output.vcf_filtered}
        """