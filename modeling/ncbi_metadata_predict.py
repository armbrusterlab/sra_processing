# turn this into a function later
# for now, use hardcoded paths and cd to the directory above the models dir

import os
import joblib
import pickle
import pandas as pd
def load_best_model(save_dir='models'):
    """
    Load the best saved model and its metadata.
    
    Parameters:
    -----------
    save_dir : str
        Directory where the model was saved
    
    Returns:
    --------
    tuple : (model, metadata)
    """
    model_path = os.path.join(save_dir, "best_model.joblib")
    metadata_path = os.path.join(save_dir, "best_model_metadata.pkl")
    
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model not found at {model_path}")
    
    # Load model
    model = joblib.load(model_path)
    
    # Load metadata
    metadata = None
    if os.path.exists(metadata_path):
        with open(metadata_path, 'rb') as f:
            metadata = pickle.load(f)
    
    return model, metadata

# model, metadata = load_best_model("models/") # in this case the best-performing model was random forest
df = pd.read_csv("/home/kcw2/data/testing/bacdive_model/1000_bacteria_metadata_transformed.tsv", sep="\t")

# vectorize the text data
from sklearn.feature_extraction.text import TfidfVectorizer

vectorizer = joblib.load('models/tfidf_vectorizer.joblib') # assumes you're currently at /home/kcw2/data/testing/bacdive_model/
selected_indices = joblib.load('models/feature_selection_indices.joblib')
y_colnames = list(joblib.load('models/y_colnames.joblib')) # keeping the list typecast to maintain compatibility with bugged outputs which produced sets instead, but it should already be a list

# Transform new data using the SAME vectorizer (NOT fit_transform!)
X = vectorizer.transform(df["joined_string"]) # X.shape returns (29, 12781)

# apply feature selection
X = X[:, selected_indices] # X.shape returns (29, 2000)

# predictions = model.predict(X)

# # the code below is robust to the possibility that the model predicts no categories/subcategories on a given joined_string
# predictions_list = [[] for i in range(len(df))]
# for i in range(len(predictions)):
#     for j in range(len(predictions[i])):
#         if predictions[i][j] == 1:
#             predictions_list[i].append(y_colnames[j])

# # separate these into category and subcategory columns

# category = []
# subcategory = []
# for t in predictions_list:
#     tt = [term.split("@@@") for term in t]
#     category.append([term_list[0] for term_list in tt])
#     subcategory.append([term_list[1] for term_list in tt])

# category_str = []
# subcategory_str = []
# for t in predictions_list:
#     tt = [term.split("@@@") for term in t]
#     category_str.append(", ".join([term_list[0] for term_list in tt]))
#     subcategory_str.append(", ".join([term_list[1] for term_list in tt]))



# # either way, add it to the dataframe.
# df["category"]=category
# df["subcategory"]=subcategory
# df["category_str"]=category_str
# df["subcategory_str"]=subcategory_str

all_models = joblib.load("models/tuning_summary.pkl")

# for model in ['random_forest', 'gradient_boosting', 'lasso', 'logistic_regression']:
for model in ['logistic_regression']: # I decided to focus on logistic regression
    print(f"Predicting with {model}...")
    predictions = all_models['models'][model].predict(X)
    predictions_list = [[] for i in range(len(df))]
    for i in range(len(predictions)):
        for j in range(len(predictions[i])):
            if predictions[i][j] == 1:
                predictions_list[i].append(y_colnames[j])

    # separate these into category and subcategory columns
    # category = []
    # subcategory = []
    # for t in predictions_list:
    #     tt = [term.split("@@@") for term in t]
    #     category.append([term_list[0] for term_list in tt])
    #     subcategory.append([term_list[1] for term_list in tt])

    category_str = []
    subcategory_str = []
    for t in predictions_list:
        tt = [term.split("@@@") for term in t]
        category_str.append(", ".join([term_list[0] for term_list in tt]))
        subcategory_str.append(", ".join([term_list[1] for term_list in tt]))
        
    df[f"category_{model}"]=category_str
    df[f"subcategory_{model}"]=subcategory_str

# finally save to file
outname = "metadata_processed.tsv"
df.to_csv(outname, sep="\t")