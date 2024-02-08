# Take a path to a directory containing fastq files and run the variant calling pipeline
input_path=$1
samples=(${input_path}"/*.fastq.gz")
reference=$2
minimap_index=$3
output_dir=$4
threads=$5

# Print usage if no arguments are given
if [ $# -eq 0 ]; then
	echo "Usage: variant_calling_pipeline.sh <input_path> <reference> <minimap_index> <output_dir> <threads>"
	exit 1
fi

mkdir -p ${output_dir}

for sample in ${samples[@]};
do
	echo "Processing sample: "$sample
	read=$sample
	sample_name=$(basename -s .fastq.gz $read)

	mkdir -p ${output_dir}/nanoplot/raw/${sample_name}_raw
	NanoPlot -t $threads --fastq $read \
		--outdir ${output_dir}/nanoplot/raw/${sample_name}_raw \
		-p ${sample_name} \
		--loglength \
		--plots dot

	mkdir -p ${output_dir}/filtlong
	filtlong --min_length 200 $read | bgzip > ${output_dir}/filtlong/${sample_name}.trimmed.fastq.gz

	read=${output_dir}/filtlong/${sample_name}.trimmed.fastq.gz

	mkdir -p ${output_dir}/nanoplot/trimmed/${sample_name}_trimmed
	NanoPlot -t $threads --fastq $read \
		--outdir ${output_dir}/nanoplot/trimmed/${sample_name}_trimmed \
		-p ${sample_name} \
		--loglength \
		--plots dot

	mkdir -p ${output_dir}/minimap2
	minimap2 -ax map-ont -t $threads $minimap_index $read | samtools sort -@ $threads -o ${output_dir}/minimap2/${sample_name}.sorted.bam
	
	samtools index ${output_dir}/minimap2/${sample_name}.sorted.bam

	bam=${output_dir}/minimap2/${sample_name}.sorted.bam
	
	sniffles --input $output_dir/minimap2/${sample_name}.sorted.bam \
			--vcf ${output_dir}/variants/SVs/${sample_name}.vcf \
			--snf ${output_dir}/variants/SVs/${sample_name}.snf \
			--reference $reference \
			--threads $threads \
			--output-rnames
	bzip2 ${output_dir}/variants/SVs/${sample_name}.snf
done
echo "Done"