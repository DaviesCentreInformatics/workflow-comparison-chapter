rule all:
	input:
		expand("results/variants/SVs/{sample}.vcf.gz", sample=config["samples"]),
		expand("results/nanoplot/raw/{sample}_raw", sample=config["samples"]),
		expand("results/nanoplot/trimmed/{sample}_trimmed", sample=config["samples"])

rule MINIMAP2:
	input:
		read = "results/filtlong/{sample}.trimmed.fastq.gz",
		index = config["minimap_index"]

	output:
		bam="results/minimap2/{sample}.sorted.bam",
		bai="results/minimap2/{sample}.sorted.bam.bai"

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
		directory("results/nanoplot/raw/{sample}_raw/")

	threads: 4

	conda:
		"environment.yaml"

	log:
		"results/logs/nanoplot/raw/{sample}.log"

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
		filtered = "results/filtlong/{sample}.trimmed.fastq.gz"

	wildcard_constraints:
		sample = "[A-Za-z0-9]+"

	threads: 4

	conda:
		"environment.yaml"

	log:
		"results/logs/filtlong/{sample}.log"

	shell: 
		"""
		filtlong --min_length 200 {input.fastq} | bgzip > results/filtlong/{wildcards.sample}.trimmed.fastq.gz
		"""

rule NANOPLOT_TRIMMED:
	input:
		fastq = "results/filtlong/{sample}.trimmed.fastq.gz"
	
	output:
		directory("results/nanoplot/trimmed/{sample}_trimmed/")

	threads: 4

	conda:
		"environment.yaml"

	log:
		"results/logs/nanoplot/trimmed/{sample}.log"

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
		bam="results/minimap2/{sample}.sorted.bam",
		bai="results/minimap2/{sample}.sorted.bam.bai",
		ref=config["reference"]

	output:
		vcf="results/variants/SVs/{sample}.vcf.gz",
		snf="results/variants/SVs/{sample}.snf.gz"

	threads: 4

	conda:
		"environment.yaml"

	log:
		"logs/variants/SVs/{sample}.log"
	
	shell:
		"""
		sniffles --input {input.bam} \
			--vcf {output.vcf} \
			--snf results/variants/SVs/{wildcards.sample}.snf \
			--reference {input.ref} \
			--threads {threads} \
			--output-rnames
		bzip2 results/variants/SVs/{wildcards.sample}.snf &> {log}
		"""