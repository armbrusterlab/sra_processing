# sra_processing
TODO: rename repository once name is finalized

# Table of contents
- [Pipeline overview](#pipeline-overview)
- [Installation](#installation)
- [Usage](#usage)
	- [Part 1: predownload](#part-1-pre-download-pipeline)
	- [Manually curate accessions](#manually-curate-accessions-to-download)
	- [Adjust environmental source predictions](#adjust-environmental-source-predictions)
 	- [Part 2: main pipeline](#part-2-main-pipeline)
- [Outputs](#outputs)
- [Acknowledgements](#acknowledgements)

# Pipeline overview
This pipeline facilitates the search for mutations of target genes within NCBI's SRA database, and summarizes the environmental sources each mutation is found in. The SRA database is broader than GenBank or RefSeq. However, SRA cannot be directly BLASTed, and the large size of SRA read datasets makes it time-intensive and compute-intensive to conduct searches across the entire database. In order to address these issues, this pipeline was designed to conduct "smart" searches limited to the most relevant datasets, without requiring genome assembly.  
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
Please modify parameters in the process block (e.g. memory, cpus) of nextflow/nextflow.config as needed. 
## Part 1: pre-download pipeline
### Params file
For variable types, please refer to the params block of predownload.nf.
#### Required parameters
* sra_query: NCBI SRA query for runs you are considering downloading. You may copy and paste the query built by the [SRA Advanced Search Builder](https://www.ncbi.nlm.nih.gov/sra/advanced/). It is completely fine to start with a broad search query, as only the metadata will be downloaded. However, you may search for specific runs or BioProjects by IDs if you already have an idea of which datasets to download.
* genome_length: For coverage calculations, the length of the genome of interest in base pairs. This may be approximate, or the average of multiple species' genome lengths. Scientific notation is supported.
* taxid: For coverage calculations, the NCBI taxonomic ID for the genus or species of interest. Only reads mapping to this taxid or its children will be counted toward coverage. Please enter this as an int, not a string, because strings will not match JSON keys.

#### Optional parameters
* coverage_threshold: In order to avoid downloading runs with very low coverage for the taxid of interest, a coverage check is conducted. Default: 30x.

### Run the pipeline
The code below will place outputs in the results/ dir. If you use the -output-dir option with a different outdir, you will need to specify the metadata/ subdir of that outdir as predownload_outputs in the params file for the main pipeline.  
```bash
cd nextflow/
nextflow run predownload.nf -params-file predownload-params.yaml
```

## Manually curate accessions to download
### Automatic suggestion of runs to download
The suggested_runs/ dir of the predownload output contains a list of the runs with the best ratio of coverage (for the user's specified taxid) to prefetch file size in megabytes for each study, ignoring any runs without predicted environmental source terms. suggested_runs/suggested_runs.txt may be used directly as run_list in the main pipeline if desired, while suggested_runs/suggested_runs_metadata.tsv can help with further narrowing down that list. This is most helpful if the user's goal is to statistically identify terms enriched in genes/mutations (which is performed at the end of the main pipeline regardless). If solely interested in discovering as many mutations as possible, it does not matter whether the runs downloaded are independent, so they don't need to all be from different studies.  
Usefully, suggested_runs/estimated_space_required.txt explains the disk space required if the user were to proceed with the suggested runs.

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

## Adjust environmental source predictions
To add field-specific custom terms or terms that the classification model missed, or to remove terms that the classification model erroneously introduced, you may use scripts/override_categorization.py.
```bash
cp "nextflow/results/metadata/environment_predictions.tsv" "nextflow/results/metadata/environment_predictions_old.tsv"

python "scripts/override_categorization.py" -p "nextflow/results/metadata/environment_predictions_old.tsv" -o "nextflow/results/metadata/environment_predictions.tsv" -a "nextflow/data/keywords_add.tsv" -m "nextflow/data/keywords_subtract.tsv"

```
Important considerations:
* You may run this script at any time, on any file with environmental source prediction columns (including files to which the predictions have been joined as metadata). HOWEVER, if running it on "environment_predictions.tsv" in between the predownload pipeline and main pipeline, please also name the new file "environment_predictions.tsv" and place it in the predownload outputs directory used for the main pipeline, at the same level as metadata_esearch.csv and metadata_pysradb.tsv.
* Be mindful of other files that may contain the out-of-date predictions.
* The -a and -m arguments are not strictly necessary, in case you only want to add terms or you only want to subtract terms.

File format:
* The -a (addition) and -m (subtraction) files are tab-separated tables with the following columns: keyword, category, subcategory
* Keywords can be written as regex patterns following the style of the Python re module. For example, word boundaries can be enforced by adding "\b" to the keyword.
* All input text is converted to lowercase, so the keywords should also be written as lowercase.
* In the -m file, generally the category and subcategory must be an exact match for existing terms in order for these terms to be removed. However, if the subcategory is provided as "*", then all terms matching the category will be removed.

## Part 2: main pipeline
### Params file
For variable types, please refer to the params block of main.nf.  
#### Required parameters
Params that should match params for the pre-download pipeline:
* genome_length
* taxid
* coverage_threshold

Other params:
* run_list: Manually curated list of runs to download, as described above.
* k2db: Path to a Kraken2 database for taxonomic classification. Please prepare this in advance. You may refer to the section on downloading a Kraken2 database.
* reference_gb: Reference sequences for target genes will be extracted from this GB file. Please refer to the "Downloading a Kraken2 database" section below for more detail.
* target_genes: space-separated list of genes in which to find variants.

#### Optional parameters
* predownload_outputs: Path to the metadata/ subdir of the predownload pipeline's outputs. By default, assumes predownload outputs are saved to the default outdir (./results/).
* read_length_threshold: Threshold in bp for whether a read is considered "short" or "long".
* fastp_additional_short and fastp_additional_long: Additional arguments for fastp and fastplong respectively.
* target_type: Indicates the field in the GB file in which to look for matching items in the target_genes list. For example, "locus_tag", "gene", or "product".
* buffer_upstream and buffer_downstream: To assist with read mapping in the breseq step, each target gene is extracted with a buffer on either side. Default for both: 900 bp.
* breseq_additional: Additional arguments for breseq.
* filter_intergenic and filter_synonymous: Filter out intergenic and/or synonymous mutations from the breseq output.
* category_colname, subcategory_colname, and terms_colname: Column names for predicted category, subcategory, and terms (where terms is a readable concatenation of the former two columns). These don't need to be changed unless you opt to use a model other than the ones provided in the models/ dir, e.g. v19 and v20. If you create your own model, refer to scripts/predict_environmental_source.py to understand how column names are assigned.
* p_adjust_method: Multiple hypothesis testing correction method; default value "fdr". Refer to [p.adjust documentation](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust) for options.
* report_all: Default value "TRUE". (Note that this is in all caps, in keeping with R boolean convention.) If not TRUE, the summary tables will only include rows with statistically significant rows.

#### Downloading a Kraken2 database
Unless you would like to build your own Kraken2 database, you may download various databases provided by Langmead et al. from [here](https://benlangmead.github.io/aws-indexes/k2). The code below produces a database at nextflow/kraken2_db/.
```bash
mkdir nextflow/kraken2_db/
curl https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz --output nextflow/kraken2_db/k2_standard_8GB.tar.gz # this is the current link for the 8 GB standard database
tar -xvzf nextflow/kraken2_db/k2_standard_8GB.tar.gz -C nextflow/kraken2_db/
rm -f nextflow/kraken2_db/k2_standard_8GB.tar.gz
```

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
cd nextflow/ # if not already at the nextflow dir
nextflow run main.nf -params-file main-params.yaml -with-report results/report.html # if using custom outdir, change report dir to match
```

# Outputs
The main output is the stats/ dir, which contains by-gene and by-mutation counts of associated environmental source terms, as well as statistical tests (2x2 Fisher's Exact Tests) to determine whether any terms are enriched for each gene/mutation. Please note that if a mutation is associated with multiple terms, it will be counted multiple times in the count tables. If you would instead like to know how many runs a particular mutation appeared in, please refer to breseq_summary_tables/mutations.tsv.  
breseq output HTMLs supporting each mutation listed in breseq_summary_tables/mutations.tsv may be found in breseq_export/. Opening this file in Excel may result in garbled text being displayed; if that is the case, please try opening it in a plaintext editor such as Notepad.  

# Acknowledgements
* Advisor: Dr. Catherine Armbruster
* Open-source projects used:
	* [pysradb](https://github.com/saketkc/pysradb)
 	* NCBI [esearch](https://eutilities.github.io/site/Quick_Start/eu_quick/#esearch) and [datasets](https://github.com/ncbi/datasets)
  	* [sra-tools](https://github.com/ncbi/sra-tools)
  	* [fastp](https://github.com/OpenGene/fastp) and [fastplong](https://github.com/OpenGene/fastplong)
  	* [kraken2](https://github.com/DerrickWood/kraken2)
  	* [breseq](https://github.com/barricklab/breseq)
