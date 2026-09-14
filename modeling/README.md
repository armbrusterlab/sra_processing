# Prepare data
## Download metadata from NCBI datasets
### Select accessions for download
For my test dataset, I wanted to randomly sample from NCBI's bacterial database.
```bash
datasets summary genome taxon bacteria --as-json-lines | dataformat tsv genome --fields accession > /home/kcw2/data/testing/bacdive_model/some_bacteria_accessions.txt # Ctrl + C to cancel partway through, leaving a file that has many, though not all, bacterial genome accessions
tail -n +2 "/home/kcw2/data/testing/bacdive_model/some_bacteria_accessions.txt" | shuf -n 1000 > /home/kcw2/data/testing/bacdive_model/1000_bacteria_accessions.txt # cut off the header and randomly pick 1000 accessions
```

### Given a text file of accessions, download metadata
```bash
accessions="/home/kcw2/data/testing/bacdive_model/1000_bacteria_accessions.txt"
raw_metadata="/home/kcw2/data/testing/bacdive_model/1000_bacteria_metadata.tsv"

datasets summary genome accession --inputfile $accessions --as-json-lines | dataformat tsv genome --fields accession,assminfo-biosample-isolation-source,assminfo-sequencing-tech,assmstats-genome-coverage,assminfo-assembly-method,assminfo-biosample-host,assminfo-biosample-host-disease,assminfo-bioproject-lineage-title,assminfo-biosample-bioproject-title,assminfo-biosample-description-title,organism-name,organism-infraspecific-strain,organism-tax-id,assmstats-contig-l50,assmstats-contig-n50,assmstats-scaffold-l50,assmstats-scaffold-n50,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > $raw_metadata
```

Not all fields listed above are necessary. The following fields are assumed to be in $raw_metadata, however:
* Assembly.BioSample.Attribute.Name
* Assembly.BioSample.Attribute.Value
* Assembly.BioProject.Lineage.Title
* Assembly.BioSample.BioProject.Title
* (possibly) Assembly.BioSample.Description.Title
* Assembly.BioSample.Isolation.source
* Assembly.BioSample.Host.disease
* Assembly.Accession

Note that spaces in column names have been replaced with periods.

### Transform the metadata
The ncbi_metadata_transform.R script joins all metadata potentially informative of environmental source into a single column, joined_string. The joined_string column is used as input for prediction.
```bash
transformed_metadata="/home/kcw2/data/testing/bacdive_model/1000_bacteria_metadata_transformed.tsv"
Rscript "scripts/modeling/ncbi_metadata_transform.R" $raw_metadata $transformed_metadata
``` 

# Predict environmental sources from data
## Train models
### Obtain training data
"scripts/modeling/isolation_sources_bacdive_2026-08-12.csv" was downloaded from [BacDive](https://bacdive.dsmz.de/isolation-sources), then processed using "scripts/modeling/wrangle_bacdive.py" to produce tidy data, located at "scripts/modeling/bacdiveReformat_2026-08-12.tsv". BacDive's level 2 tags were concatenated with their corresponding level 1 tags to produce unique labels upon which a multi-label classification model could be built. These concatenated labels may be found in the joined_1_2 column of "scripts/modeling/bacdiveReformat_2026-08-12.tsv". To reconstitute a list of concatenated level 1 and level 2 tags, split this column by the "###" delimiter. To separate level 1 and level 2 tags, the strings may be further split by the "@@@" delimiter.  
Note that this data is not sourced from NCBI like the test data is.

### Clean training data
In "scripts/modeling/bacdive_trainModel.py", various stopwords are provided, including names of bacteria and of countries, to eliminate false associations between these names and environmental sources. The input text was also cleaned of additional characters that could interfere with string parsing.  
I made the decision to remove certain BacDive labels from the input data entirely. I removed everything under the Host level 1 category because it was ambiguous (for example, in some cases "dog feces" would be tagged with both "host, mammal" and "host, microbial"). I also removed the "infection, patient" category because it was strongly associated with bodily products that don't necessarily indicate that the host was infected, for example "urine".

### Build models
For the predictor data, the "Isolation source" column of "scripts/modeling/bacdiveReformat_2026-08-12.tsv" was processed with a TF-IDF text vectorizer. As for the response data, the joined_1_2 column of "scripts/modeling/bacdiveReformat_2026-08-12.tsv" was multi-hot encoded. I used an 80:20 train-test split.  
Feature selection was performed using a hybrid of chi-square and random forest to select the 2000 most informative words among the predictors.  
Optuna was used for 5-fold cross-validation of various models using the feature-selected data. I built and tested the following kinds of classification models: random forest, gradient boosting, Lasso regression, and more generally logistic regression.
```bash
model_dir="vx" # e.g. v1, v2... the model-building script is currently hard-coded and will produce outputs in the current directory
mkdir $model_dir
cd $model_dir

python "~/metadata-magnet/scripts/modeling/bacdive_trainModel.py"
```

## Test models
 I tested the models on $transformed_metadata.
 ```bash 
 # still in $model_dir
 python "~/metadata-magnet/scripts/modeling/ncbi_metadata_predict.py"
 ```
 In doing so, I found that over many iterations of the models, the logistic regression model tends to perform the best, so I decided to focus on that. I tested combinations of metrics (Jaccard vs precision; macro vs weighted averaging strategy) and settled on two models: Jaccard macro for a less conservative model that will occasionally assign environmental sources that don't apply, and precision macro for a more conservative model that is usually accurate with its labeling but will sometimes fail to classify isolation sources that Jaccard macro would be able to. 

 # Future directions
 The models could be improved by adding to the training data based on shortcomings of the prediction upon $transformed_metadata. For example, if the word "bed" does not appear enough in the training data to be associated with the "built environment" subcategory, then adding some lines with the isolation source and user-verified labels would reinforce that connection in the model. This would require extensive manual review of the predictions, but it may be worth investing the time to do so. Additionally, many entries in the training data had isolation sources but no tags and so I had to discard them; it is possible that more tags will be added later, as the BacDive team tags these entries manually.
 I considered using n-grams rather than single words when vectorizing the text. However, the ordering of words may be somewhat arbitrary because I join data from unrelated columns into a single column to produce the training data.