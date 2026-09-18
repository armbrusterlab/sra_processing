#!/usr/bin/env python

import pandas as pd
import re
import joblib
import argparse

def predict_env_source(metadata_file, model_dir, outname="env_source_predictions.tsv"):
    print("Running join_strings...")
    df = join_strings(metadata_file)
    print("Running predict_on_metadata...")
    predict_on_metadata(df, model_dir, outname)


def join_strings(f):
    df = pd.read_csv(f, sep="\t")

    # there's some variation in column names, but columns with names containing any of these keywords should be relevant
    keywords=["study_title", "isolation", "environment", "organism part", "tissue", "env_biome", "disease"]
    matching_cols = [col for col in df.columns if any(substring in col for substring in keywords)]

    # for study_title column specifically, only extract the part following "from" so as to avoid using misleading strings from the full title
    if "study_title" in matching_cols:
        df["study_title_from"] = df["study_title"].str.extract(r"from (.*)", flags=re.IGNORECASE)
        idx = matching_cols.index("study_title")
        matching_cols[idx]="study_title_from"

    new_df = pd.DataFrame()
    new_df['run_accession'] = df['run_accession'] # df is from pysradb, which has the run_accession column

    new_df['text_for_prediction'] = (
        df[matching_cols]
        .fillna('') # NA's become empty strings
        .agg(' '.join, axis=1)  # join with space
    )

    # df.insert(1, "text_for_prediction", df.pop("text_for_prediction"))

    return new_df


def predict_on_metadata(df, model_dir, outname):
    vectorizer = joblib.load(f'{model_dir}/tfidf_vectorizer.joblib') # assumes you're currently at /home/kcw2/data/testing/bacdive_model/
    selected_indices = joblib.load(f'{model_dir}/feature_selection_indices.joblib')
    y_colnames = joblib.load(f'{model_dir}/y_colnames.joblib')

    # transform data using the same vectorizer as was used to build the model
    X = vectorizer.transform(df["text_for_prediction"])

    # apply feature selection
    X = X[:, selected_indices] # X.shape returns (29, 2000)

    all_models = joblib.load(f"{model_dir}/tuning_summary.pkl")

    for model in all_models['models'].keys():
        print(f"Predicting with {model}...")
        predictions = all_models['models'][model].predict(X)
        predictions_list = [[] for i in range(len(df))]
        for i in range(len(predictions)):
            for j in range(len(predictions[i])):
                if predictions[i][j] == 1:
                    predictions_list[i].append(y_colnames[j])

        # separate these into category and subcategory columns
        category = []
        subcategory = []
        terms = []
        for t in predictions_list:
            tt = [term.split("@@@") for term in t]
            category.append([term_list[0] for term_list in tt])
            subcategory.append([term_list[1] for term_list in tt])
            terms.append(", ".join([term.replace("@@@", ": ") for term in t]))

        df[f"category_{model}"]=category
        df[f"subcategory_{model}"]=subcategory
        df[f"terms_{model}"]=terms

        # # the code below saves the predictions as strings rather than lists
        # category_str = []
        # subcategory_str = []
        # for t in predictions_list:
        #     tt = [term.split("@@@") for term in t]
        #     category_str.append(", ".join([term_list[0] for term_list in tt]))
        #     subcategory_str.append(", ".join([term_list[1] for term_list in tt]))
            
        # df[f"category_{model}"]=category_str
        # df[f"subcategory_{model}"]=subcategory_str

        df.insert(2, f"terms_{model}", df.pop(f"terms_{model}"))
        df.insert(3, f"category_{model}", df.pop(f"category_{model}"))
        df.insert(4, f"subcategory_{model}", df.pop(f"subcategory_{model}"))


    df.to_csv(outname, sep="\t", index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to predict environmental sources from a metadata file.")

    parser.add_argument("-m", "--metadata_file", type=str, help='Metadata file, e.g. from pysradb.')
    parser.add_argument("-d", "--model_dir", type=str, help="Directory containing files for classification model.")
    parser.add_argument("-o", "--outname", type=str, default="env_source_predictions.tsv", help="Output filename.") # optional

    args = parser.parse_args()
    print(args.metadata_file, args.model_dir, args.outname)
    predict_env_source(args.metadata_file, args.model_dir, args.outname)