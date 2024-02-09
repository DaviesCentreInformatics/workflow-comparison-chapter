rule all:
	input:
		expand("snakemake_results/variants/SVs/{sample}.vcf.gz", sample=config["samples"]),
		expand("snakemake_results/nanoplot/raw/{sample}_raw", sample=config["samples"]),
		expand("snakemake_results/nanoplot/trimmed/{sample}_trimmed", sample=config["samples"])

rule MINIMAP2:
	input:
		read = "snakemake_results/filtlong/{sample}.trimmed.fastq.gz",
		index = config["minimap_index"]

	output:
		bam="snakemake_results/minimap2/{sample}.sorted.bam",
		bai="snakemake_results/minimap2/{sample}.sorted.bam.bai"

	threads: 8
	
	conda:
		"environment.yaml"

	log:
		"logs/minimap2/{sample}.log"
	
	shell:
		"""
		minimap2 -ax map-ont -t {threads} {input.index} {input.read} | samtools sort -@ {threads} -o {output.bam} &> {log}

		samtools index {output.bam}
		"""

rule NANOPLOT_RAW:
	input:
		fastq = "input_data/raw_fastq/{sample}.fastq.gz"
	
	output:
		directory("snakemake_results/nanoplot/raw/{sample}_raw/")

	threads: 4

	conda:
		"environment.yaml"

	log:
		"snakemake_results/logs/nanoplot/raw/{sample}.log"

	shell: 
		"""
		outprefix=$(basename -s .fastq.gz {input.fastq})_raw
		NanoPlot -t {threads} \
			--fastq {input.fastq} \
			--outdir {output} \
			-p $outprefix \
			--loglength \
			--plots dot &> {log}
		"""

rule FILTLONG:
	input:
		fastq = "input_data/raw_fastq/{sample}.fastq.gz"
	
	output:
		filtered = "snakemake_results/filtlong/{sample}.trimmed.fastq.gz"

	wildcard_constraints:
		sample = "[A-Za-z0-9]+"

	threads: 4

	conda:
		"environment.yaml"

	log:
		"snakemake_results/logs/filtlong/{sample}.log"

	shell: 
		"""
		filtlong --min_length 200 {input.fastq} | bgzip > snakemake_results/filtlong/{wildcards.sample}.trimmed.fastq.gz
		"""

rule NANOPLOT_TRIMMED:
	input:
		fastq = "snakemake_results/filtlong/{sample}.trimmed.fastq.gz"
	
	output:
		directory("snakemake_results/nanoplot/trimmed/{sample}_trimmed/")

	threads: 4

	conda:
		"environment.yaml"

	log:
		"snakemake_results/logs/nanoplot/trimmed/{sample}.log"

	shell: 
		"""
		outprefix=$(basename -s .fastq.gz {input.fastq})_trimmed
		NanoPlot -t {threads} \
			--fastq {input.fastq} \
			--outdir {output} \
			-p $outprefix \
			--loglength \
			--plots dot &> {log}
		"""

rule SNIFFLES2:
	input:
		bam="snakemake_results/minimap2/{sample}.sorted.bam",
		bai="snakemake_results/minimap2/{sample}.sorted.bam.bai",
		ref=config["reference"]

	output:
		vcf="snakemake_results/variants/SVs/{sample}.vcf.gz",
		snf="snakemake_results/variants/SVs/{sample}.snf.bz2"

	threads: 4

	conda:
		"environment.yaml"

	log:
		"logs/variants/SVs/{sample}.log"
	
	shell:
		"""
		sniffles --input {input.bam} \
			--vcf {output.vcf} \
			--snf snakemake_results/variants/SVs/{wildcards.sample}.snf \
			--reference {input.ref} \
			--threads {threads} \
			--output-rnames
		bzip2 snakemake_results/variants/SVs/{wildcards.sample}.snf &> {log}
		"""