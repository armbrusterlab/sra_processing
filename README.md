# sra_processing
TODO: rename repository once name is finalized

# Table of contents
- [Pipeline overview](#pipeline-overview)
- [Installation](#installation)
- [Usage](#usage)
	- [Part 1: predownload](#part-1-pre-download-pipeline)
	- [Manually curate accessions](#manually-curate-accessions-to-download)
 	- [Part 2: main pipeline](#part-2-main-pipeline)
- [Outputs](#outputs)
- [Acknowledgements](#acknowledgements)

# Pipeline overview
This pipeline facilitates the search for mutations of target genes within NCBI's SRA database. The SRA database is broader than GenBank or RefSeq. However, SRA cannot be directly BLASTed, and the large size of SRA read datasets makes it time-intensive and compute-intensive to conduct searches across the entire database. In order to address these issues, this pipeline was designed to conduct "smart" searches limited to the most relevant datasets, without requiring genome assembly.  
<img width="500" alt="image" src="https://github.com/user-attachments/assets/76e9b108-eb43-4ed5-943e-453522ce4ab3" />


# Installation
Developed in Ubuntu 24.04.2 LTS (GNU/Linux 6.8.0-44-generic x86_64) for Nextflow v25.10.4. If you encounter any issues running the pipeline, please ensure that your Nextflow version matches this version, as Nextflow is under active development and other versions may not be compatible with this pipeline.
```bash
conda create --name nf-env bioconda::nextflow==25.10.4
conda activate nf-env
export NXF_SYNTAX_PARSER=v2

nextflow -v # confirm that the version being used is 25.10.4
```
Any other dependencies are handled by Nextflow. Refer to envs/envs.yml.

# Usage
## Part 1: pre-download pipeline
### Params file
For variable types, please refer to the params block of predownload.nf.
* sra_query: NCBI SRA query for runs you are considering downloading. You may copy and paste the query built by the [SRA Advanced Search Builder](https://www.ncbi.nlm.nih.gov/sra/advanced/). It is completely fine to start with a broad search query, as only the metadata will be downloaded. However, you may search for specific runs or BioProjects by IDs if you already have an idea of which datasets to download.
* genome_length: For coverage calculations, the length of the genome of interest in base pairs. This may be approximate, or the average of multiple species' genome lengths. Scientific notation is supported.
* taxid: For coverage calculations, the NCBI taxonomic ID for the genus or species of interest. Only reads mapping to this taxid or its children will be counted toward coverage. Please enter this as an int, not a string, because strings will not match JSON keys.
* coverage_threshold: In order to avoid downloading runs with very low coverage for the taxid of interest, a coverage check is conducted. Default: 30x.

### Run the pipeline
The code below will place outputs in the results/ dir. If you use the -output-dir option with a different outdir, you will need to specify the metadata/ subdir of that outdir as predownload_outputs in the params file for the main pipeline.  
```bash
nextflow run predownload.nf -params-file predownload-params.yaml
```

## Manually curate accessions to download
### Disk space considerations
SRA runs are very large. Please refer to metadata in metadata/coverage_check/ from the pre-download pipeline in order to select runs to download in the main pipeline. To estimate the disk space required by a set of runs, sum their size_megabytes values from the metadata and multiply that by 25 or 30. (The actual fastq files are about [7 times the size of the prefetch files](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump/633360aa302b9f2b6e8ceda7c99dde07f4f20e2e#extract-fastq-files-from-the-sra---accession), and due to [Nextflow's limitations](https://github.com/nextflow-io/nextflow/issues/452), the 4 processes' worth of redundant files are not deleted once no longer needed by the pipeline. To save space, please use `nextflow clean` to remove intermediate files.)  
The conda environment used by the Nextflow pipeline requires about 2.8 GB.

### Required file format
The curated list of accessions must be a single-column text file with one run accession on each line and no header. As an example of the correct formatting, refer to results_example/metadata/coverage_check/coverage_pass.txt.  
In curating the accessions, you may remove accessions that appeared in the initial query, but please do not add accessions as they will lack metadata used downstream in the pipeline.  
If there are Windows-style return characters in run_list, it may cause issues with parsing. In that case, please sanitize the file as demonstrated below.
```bash
run_list="my_run_list.txt"
sed -i 's/\r$//' $run_list
```

## Part 2: main pipeline
### Params file
For variable types, please refer to the params block of main.nf.  
Params that should match params for the pre-download pipeline:
* genome_length
* taxid
* coverage_threshold

Other params:
* run_list: Manually curated list of runs to download, as described above.
* predownload_outputs: Path to the metadata/ subdir of the predownload pipeline's outputs. By default, assumes predownload outputs are saved to the default outdir (./results/).
* read_length_threshold: Threshold in bp for whether a read is considered "short" or "long".
* fastp_additional_short and fastp_additional_long: Additional arguments for fastp and fastplong respectively.
* k2db: Path to a Kraken2 database for taxonomic classification. Please prepare this in advance.
* reference_gb: Reference sequences for target genes will be extracted from this GB file. Please refer to the section below for more detail.
* target_genes: space-separated list of genes in which to find variants.
* target_type: Indicates the field in the GB file in which to look for matching items in the target_genes list. For example, "locus_tag", "gene", or "product".
* buffer_upstream and buffer_downstream: To assist with read mapping in the breseq step, each target gene is extracted with a buffer on either side. Default for both: 900 bp.
* breseq_additional: Additional arguments for breseq.

#### Obtaining a GB file
reference_gb may either be downloaded from a website or obtained from NCBI datasets. See below for the latter approach.
```bash
datadir="." # in this case, the repo dir
acc="GCF_000006765.1"

datasets download genome accession $acc --include gbff --filename ${datadir}/ncbi_download.zip
unzip -d $datadir ${datadir}/ncbi_download.zip
rm -f ${datadir}/ncbi_download.zip
						
find ${datadir}/ncbi_dataset -type f -name genomic.gbff # prints the filename; paste this into the params file as reference_gb
```

### Run the pipeline
The outdir does not need to match that of the predownload pipeline.
```bash
nextflow run main.nf -params-file main-params.yaml -with-report results/report.html
```

# Outputs
The main output is breseq_summary_tables/breseq_summary_withMetadata.tsv, which joins metadata from the predownload pipeline to the mutations found by breseq relative to target gene sequences in reference_gb. (The missing coverage and new junction outputs are also aggregated, but do not include metadata.) breseq output HTMLs supporting each mutation may be found in breseq_export/.

# Acknowledgements
TODO
