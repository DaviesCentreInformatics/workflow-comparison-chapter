#!/usr/bin/env nextflow
// Set the DSL syntax to version 2
nextflow.enable.dsl=2

/*
 * PARAMETERS
 */
params.samplesheet = null
params.reference = null
params.reference_index = null
params.minimap_index = null
params.outdir = null


/*
 * CHECK INPUTS
 */

if (params.samplesheet == null) {
	error "Please provide a sample sheet using the `--samplesheet` flag"
	System.exit(1)
}
if (params.minimap_index == null) {
	error "Please provide a minimap_index using the `--minimap_index` flag"
	System.exit(1)
}

if (params.reference == null) {
	error "Please provide a reference genome using the `--reference` flag"
	System.exit(1)
}

if (params.reference_index == null) {
	error "Please provide a reference genome index using the `--reference_index` flag"
	System.exit(1)
}

if (params.minimap_index == null) {
	error "Please provide a minimap index using the `--minimap_index` flag"
	System.exit(1)
}
if (params.outdir == null) {
	error "Please provide an output directory using the `--outdir` flag"
	System.exit(1)
}

log.info """\
N E X T F L O W  --  V A R I A N T   C A L L I N G  --  O N T
=============================================================

Pipeline Parameters:
		Sample sheet:     ${params.samplesheet}
		Output directory: ${params.outdir}
		Reference genome: ${params.reference}
		Reference index:  ${params.reference_index}
		Minimap index:    ${params.minimap_index}

=============================================================
"""

// Create input channel from the sample sheet
Channel.fromPath(params.samplesheet, checkIfExists: true)
					 .splitCsv(header: true)
		             .map ( row -> tuple(row.sampleID, row.fastq) )
					 .set { samples }

process NANOPLOT_RAW {
	tag "$sampleID"
	label "process_low", "error_retry"

	publishDir "$params.outdir/nanoplot", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("trimmed") > 0 ) "trimmed/$filename"
			else "raw/$filename" }

	input:
	tuple val(sampleID), path(fastq)

	output:
	path "${sampleID}_raw/*", emit: report

	script:
	outdir = "${sampleID}_raw"
	outprefix = "${sampleID}_raw"

	"""
	NanoPlot -t ${task.cpus} \
		--fastq ${fastq} \
		--outdir $outdir \
		-p $outprefix \
		--loglength \
		--plots dot
	"""
}

process FILTLONG {
	tag "$sampleID"
	label 'process_low', 'error_retry'

	publishDir "$params.outdir/filtlong", mode: 'copy'

	input:
	tuple val(sampleID), path(fastq)

	output:
	tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz") , emit: result_tuple
	
	script:
	"""
	filtlong --min_length 200 ${fastq} | bgzip > ${sampleID}.trimmed.fastq.gz
	"""
}

process NANOPLOT_FILTERED {
	tag "$sampleID"
	label "process_low", "error_retry"

	publishDir "$params.outdir/nanoplot/trimmed", mode: 'copy'

	input:
	tuple val(sampleID), path(fastq)

	output:
	path "${sampleID}_trimmed/*", emit: report

	script:
	outdir = "${sampleID}_trimmed"
	outprefix = "${sampleID}_trimmed_"

	"""
	NanoPlot -t ${task.cpus} \
		--fastq ${fastq} \
		--outdir $outdir \
		-p $outprefix \
		--loglength \
		--plots dot
	"""
}

process MINIMAP2 {
	tag "$sampleID"
	label "process_high"

	publishDir "$params.outdir/minimap2", mode: 'copy'

	input:
	tuple val( sampleID ), path( fastq )
	path reference

	output:
	tuple val( sampleID ), path( "${sampleID}.sorted.bam" ), path( "${sampleID}.sorted.bam.bai" ), emit: mapped_tuple

	script:
	"""
	minimap2 -ax map-ont --MD -t $task.cpus $reference $fastq | samtools sort -@ $task.cpus -o ${sampleID}.sorted.bam

	samtools index ${sampleID}.sorted.bam
	"""
}

process SPLITBAM {
	tag "$sampleID"
	label "process_low"

	//publishDir "$params.outdir/variants/splitBams", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)

	output:
    tuple val(sampleID), path("${sampleID}.*.bam"), path("${sampleID}.*.bam.bai"), emit: split_bam

	shell:
    '''
    samtools idxstats !{bam} | cut -f 1 | grep -v '*' > !{sampleID}.chromosomes.txt
    while IFS= read -r line; do
        samtools view -b !{bam} ${line} > !{sampleID}.${line}.bam ;
        samtools index !{sampleID}.${line}.bam
    done < !{sampleID}.chromosomes.txt
    '''
}

process SNIFFLES2 {
	tag "$sampleID"
	label "process_medium", "error_retry"

	publishDir "$params.outdir/variants/SVs", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)
	path fa
	path faidx

	output:
	path "*"
	// Need to create a more specific output definition

	// TODO: Make this pattern more generalisable or configurable.
	script:
	"""
	sniffles --input $bam \
		--vcf ${sampleID}.vcf.gz \
		--snf ${sampleID}.snf \
		--reference $fa \
		--threads ${task.cpus} \
		--output-rnames
	bzip2 ${sampleID}.snf
	"""
}


workflow {
	NANOPLOT_RAW(samples)
	FILTLONG(samples)
	filtered_reads = FILTLONG.out.result_tuple
	NANOPLOT_FILTERED(filtered_reads)

	MINIMAP2(filtered_reads, params.minimap_index)
	mapped = MINIMAP2.out.mapped_tuple
	
	SNIFFLES2(mapped, params.reference, params.reference_index)
}