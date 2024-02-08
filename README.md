# workflow-comparison-chapter

## Nextflow

Files necessary for the Nextflow version of the pipeline:

- `conf/base.config`
- `nextflow-dockerexample.config`
- `main.nf`
- `nextflow.config`


Run the Nextflow pipeline with the following command:

```bash
conda activate nextflow
nextflow run 
```

## Snakemake

Files necessary for the Snakemake version of the pipeline:

- `conf/snakemake-config`
- `environment.yaml`
- `Snakefile`


Run the Snakemake pipeline with the following command:

```bash
conda activate snakemake
snakemake -c 8 --use-conda --configfile conf/snakemake-config.yaml
```

## Shell

File necessary for the Shell version of the pipeline:

- `variant_calling_pipeline.sh`

__Note:__ If you do not have the required software installed you can use 
`environment.yaml` to create a conda environment with all the necessary software.

If you already have the necessary software, you can skip this step.

```bash
conda env create -f environment.yaml
```

Run the Shell pipeline with the following command:

```bash
conda activate variant # or whatever you named the environment
bash variant_calling_pipeline.sh ./input_data/raw_fastq ./input_data/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.fa ./input_data/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.mmi shell_results 8
```