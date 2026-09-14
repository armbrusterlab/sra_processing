import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
# from sklearn.datasets import make_multilabel_classification
from sklearn.multioutput import MultiOutputClassifier
from sklearn.linear_model import LogisticRegression

myfile = r"C:\Users\achro\OneDrive\Desktop\CMU\Spring 2025\Armbruster Lab research\bacdiveReformat_2026-08-12.tsv"
df = pd.read_csv(myfile, sep="\t")
print(len(df)) # 63296

# start filtering
df = df[df["Isolation source"].isna() == False]
print(len(df)) # 58625

df = df[df["joined_1_2"].isna() == False]
print(len(df)) # 43314

# reconvert from string to list
df["joined_1_2"] = [
    row["joined_1_2"].split("###")
    for _, row in df.iterrows()
]

terms = {} # frequency dictionary
for t in df["joined_1_2"]:
    for term in t:
        terms[term] = terms.get(term, 0) + 1

freqs = list(terms.values())

plt.hist(freqs, bins='auto')
plt.xlabel("Count")
plt.ylabel("Frequency")
plt.title("Histogram of counts")
# plt.show()

# histogram x axis range is wide so I want to zoom in on the region that shows less-frequent terms
terms_infrequent = {k: v for k, v in terms.items() if v < 2000}
plt.hist(list(terms_infrequent.values()), bins='auto')
plt.xlabel("Count")
plt.ylabel("Frequency")
plt.title("Histogram of counts")
# plt.show()

# large spike somewhere below count=250 so zoom in even further
terms_infrequent_2 = {k: v for k, v in terms_infrequent.items() if v < 250}
plt.hist(list(terms_infrequent_2.values()), bins='auto')
plt.xlabel("Count")
plt.ylabel("Frequency")
plt.title("Histogram of counts")
# plt.show()

# zoom in even further
terms_infrequent_3 = {k: v for k, v in terms_infrequent_2.items() if v < 50}
plt.hist(list(terms_infrequent_3.values()), bins='auto')
plt.xlabel("Count")
plt.ylabel("Frequency")
plt.title("Histogram of counts")
# plt.show()

# zoom in even further...
terms_infrequent_4 = {k: v for k, v in terms_infrequent_3.items() if v < 10}
plt.hist(list(terms_infrequent_4.values()), bins='auto')
plt.xlabel("Count")
plt.ylabel("Frequency")
plt.title("Histogram of counts")
# plt.show()
# there are 30 terms that appear only once. Absolutely can't train anything on these labels, so we need to throw them out.

# let's filter out the infrequent terms
len(terms) # 117
len({k:v for k, v in terms.items() if v >= 10}) #69
len({k:v for k, v in terms.items() if v >= 20}) # 65; don't lose much by requiring 20 rather than 10
len({k:v for k, v in terms.items() if v >= 50}) # 54

# I've read that for good classification you want at least 50 items per label, but it might still be possible with 20
# Even if there are 65 different labels, training the models should be reasonable
minSize = 20
terms_blacklist = set([k for k in terms.keys() if terms[k] < minSize])
terms_keep = sorted([k for k in terms.keys() if terms[k] >= minSize])

# filter out infrequent terms from the joined_1_2 column
df["joined_1_2"] = [
    list(set(row["joined_1_2"]) - terms_blacklist)
    for _, row in df.iterrows()
]

# some rows may now have empty lists, so remove these rows
df = df[df["joined_1_2"].str.len() > 0] # even though the items are lists, pandas's str.len() function can get list lengths
print(len(df)) # 43283; compared to 43314, didn't remove too much

# now convert the response column, joined_1_2, into binary yes/no columns, one per label
y = pd.DataFrame()

for label in terms_keep:
    print(f"Processing {label}...")
    y[label] = [1 if label in t else 0 for t in df["joined_1_2"]]

# now vectorize the isolation source text data to make the feature columns, and perform feature selection
# https://www.geeksforgeeks.org/nlp/text-classification-using-scikit-learn-in-nlp/
# Initialize TF-IDF Vectorizer
vectorizer = TfidfVectorizer(stop_words='english', max_df=0.7) # first-pass "feature selection"; remove words that appear in >70% of the records since they're not helpful
# I could use min_df but opted not to, as many groups are very small

# Transform the text data to feature vectors
# first, convert the ### delimiter in the isolation source strings into spaces so that the words on either side of the delimiter will be processed properly
# also remove the <I> and </I> for italicization
iso_source = [
    row["Isolation source"].replace("###", " ").replace("<I>", "").replace("</I>", "")
    for _, row in df.iterrows()
]
# X = vectorizer.fit_transform(df['Isolation source'])
X = vectorizer.fit_transform(iso_source)

# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# wait I want to know the indices corresponding to the train-test split, so I need to split the indices instead of the df directly
indices = np.arange(len(df))

train_idx, test_idx = train_test_split(
    indices, test_size=0.2, random_state=42
)

X_train = X[train_idx]
X_test  = X[test_idx]
y_train = y.iloc[train_idx]
y_test  = y.iloc[test_idx]


clf = MultiOutputClassifier(LogisticRegression()).fit(X_train, y_train) # this is an example; you can put any type of model into the MultiOutputClassifier
# Logistic regression is pretty simple, so it takes a few seconds (~10 seconds maybe?) to finish training


# for i in range(20):
#     print(df["Isolation source"].iloc[i], df["joined_1_2"].iloc[i])
#     for a, b in zip(list(clf.predict(X_test[i])[0]), terms_keep):
#         if a == 1:
#             print(b)
#     print()

def jaccard_similarity(list1, list2):
    intersection = len(set(list1).intersection(set(list2)))
    union = len(set(list1).union(set(list2)))
    return intersection / union if union != 0 else 0.0

mismatch_count = 0
jaccard_sum = 0
n_to_test = 100

# I've discovered that if I copy and paste code into the terminal, I can't have empty lines inside a loop, otherwise it thinks I ended the loop and am entering unrelated indented lines after the empty space
for i in range(n_to_test):
    orig_idx = test_idx[i]                 # original df index
    # print(df["Isolation source"].iloc[orig_idx])
    # print(df["joined_1_2"].iloc[orig_idx])
    preds = clf.predict(X_test[i])[0]      # prediction for this training sample
    # print("Predictions:")
    assignments = []
    for a, b in zip(list(preds), terms_keep):
        if a == 1:
            assignments.append(b)
    matched = (sorted(df["joined_1_2"].iloc[orig_idx]) == sorted(assignments))
    # print(f"Perfect match? {matched}")
    if not matched:
        print("Mismatch")
        print(df["Isolation source"].iloc[orig_idx])
        print(iso_source[orig_idx])
        print(sorted(df["joined_1_2"].iloc[orig_idx]))
        print(sorted(assignments))
        similarity = jaccard_similarity(df["joined_1_2"].iloc[orig_idx], assignments)
        print(f"Jaccard Similarity: {similarity:.2f}")
        jaccard_sum += similarity
        print()
        mismatch_count += 1
    # print()

print(f"{mismatch_count} / {n_to_test} weren't perfect matches") # 25/100 not perfect matches; greatly improved from 50/100 from before I cleaned up the input text
print(f"Average Jaccard similarity among mismatches: {jaccard_sum/mismatch_count}")
print(f"Average Jaccard similarity if including perfect matches: {((n_to_test - mismatch_count) + jaccard_sum) / n_to_test}") # this is the metric I can use to evaluate model performance
# the way to handle this in the real code: calculate Jaccard for each one, and it'll be 1 if there's a perfect match.
# Increment the Jaccard sum with each Jaccard, and at the end, divide by the total size of the test (or validation I guess) dataset.

# >>> print(f"{mismatch_count} / {n_to_test} weren't perfect matches") # 25/100 not perfect matches; greatly improved from 50/100 from before I cleaned up the input text
# 25 / 100 weren't perfect matches
# >>> print(f"Average Jaccard similarity among mismatches: {jaccard_sum/mismatch_count}")
# Average Jaccard similarity among mismatches: 0.41
# >>> print(f"Average Jaccard similarity if including perfect matches: {((n_to_test - mismatch_count) + jaccard_sum) / n_to_test}")
# Average Jaccard similarity if including perfect matches: 0.8525

# look at the features (vectorized text data)
print(f"Dimensions: {X.shape}") # (43283, 12781); there are over 12k features to handle
# never mind on the below; I was going to visualize a heatmap, but this might not be feasible with how many features there are
X_df = pd.DataFrame(X.todense()) # X itself is a CSR matrix
# correlation_matrix = X_df.corr()

# I'll wrap the dimensionality reduction stuff in a function so that I can collapse it and not run it again

def dimensionality_reduction():
    # dimensionality reduction of features: PCA, t-SNE
    # this can help with looking for outliers
    import seaborn as sns
    from sklearn.decomposition import PCA

    # we already have the PCA features (X_df), just need to make the PCA target
    # I want to color by the level 1 categories
    pca_target = [
        "###".join(sorted(list(set(row["Category 1"].split("###")))))
        for _, row in df.iterrows()
    ]

    # # 1. Choose the model class
    # model = PCA(n_components=2)  # 2 components

    # # 2. Instantiate the model and fit it to the data
    # model.fit(X_df)

    # # 3. Transform the data to two dimensions
    # X_2D = model.transform(X_df)

    # # 4. Add PCA results to the DataFrame
    # df['PCA1'] = X_2D[:, 0]
    # df['PCA2'] = X_2D[:, 1]

    # # # 5. Plot the PCA results
    # df["pca_target"] = pca_target
    # sns.lmplot(x="PCA1", y="PCA2", hue='pca_target', data=df, fit_reg=False)
    # plt.show()
    # # lmplot works since 'pca_target' is categorical
    # # no clear clusters

    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler

    # Standardize features (important for t-SNE)
    # X_scaled = StandardScaler().fit_transform(X_df)

    # wait, StandardScaler does work on CSR matrix, but need with_mean = False to keep the data sparse and avoid memory issues
    scaler = StandardScaler(with_mean=False)
    X_scaled = scaler.fit_transform(X)

    # Apply t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    X_tsne = tsne.fit_transform(X_scaled) # took a few minutes to finish running

    # Plotting the results
    # codes, uniques = pd.factorize(df["pca_target"])

    # plt.figure(figsize=(8, 6))
    # plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=codes, cmap='viridis')
    # plt.colorbar()
    # plt.title('t-SNE')
    # plt.xlabel('t-SNE Component 1')
    # plt.ylabel('t-SNE Component 2')
    # plt.show()

    # # for discrete colors
    # from matplotlib.colors import ListedColormap

    # codes, uniques = pd.factorize(df["pca_target"])
    # cmap = ListedColormap(plt.cm.tab20(np.linspace(0, 1, len(uniques))))

    # plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=codes, cmap=cmap, s=100)
    # plt.colorbar()

    # # with legend instead of colorbar
    # codes, uniques = pd.factorize(df["pca_target"])

    # from matplotlib.colors import ListedColormap
    # cmap = ListedColormap(plt.cm.tab20(np.linspace(0, 1, len(uniques))))

    # plt.figure(figsize=(8, 6))
    # plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=codes, cmap=cmap, s=1)
    # plt.title('t-SNE')
    # plt.xlabel('t-SNE Component 1')
    # plt.ylabel('t-SNE Component 2')

    # from matplotlib.lines import Line2D

    # handles = [
    #     Line2D(
    #         [], [], 
    #         marker='o', 
    #         color=cmap(i), 
    #         linestyle='',
    #         markersize=10,
    #         label=uniques[i]
    #     )
    #     for i in range(len(uniques))
    # ]

    # plt.legend(handles=handles, title="Categories", bbox_to_anchor=(1.05, 1), loc='upper left')
    # plt.tight_layout()
    # plt.show()

    # # tab20 doesn't have enough colors, so go back to viridis
    # codes, uniques = pd.factorize(pca_target)
    # n_classes = len(uniques)

    # viridis = plt.cm.viridis
    # colors = viridis(np.linspace(0, 1, n_classes))

    # plt.figure(figsize=(8, 6))
    # plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=codes, cmap=viridis, s=1)
    # plt.title('t-SNE')
    # plt.xlabel('t-SNE Component 1')
    # plt.ylabel('t-SNE Component 2')

    # from matplotlib.lines import Line2D

    # handles = [
    #     Line2D(
    #         [], [], 
    #         marker='o', 
    #         color=colors[i], 
    #         linestyle='',
    #         markersize=10,
    #         label=uniques[i]
    #     )
    #     for i in range(n_classes)
    # ]

    # plt.legend(handles=handles, title="Categories", bbox_to_anchor=(1.05, 1), loc='upper left')
    # plt.tight_layout()
    # plt.show()


    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import ListedColormap

    # number of categories
    codes, uniques = pd.factorize(df["pca_target"])
    n_classes = len(uniques)

    # build a large discrete palette
    colors = (
        plt.cm.tab20(np.linspace(0, 1, 20)).tolist() +
        plt.cm.tab20b(np.linspace(0, 1, 20)).tolist() +
        plt.cm.tab20c(np.linspace(0, 1, 20)).tolist() +
        sns.color_palette("husl", max(0, n_classes - 60))  # fill remaining colors
    )

    # truncate or extend to exactly n_classes
    colors = colors[:n_classes]

    # convert to ListedColormap for scatter
    cmap = ListedColormap(colors)

    plt.figure(figsize=(8, 6))
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=codes, cmap=cmap, s=1)
    plt.title('t-SNE')
    plt.xlabel('t-SNE Component 1')
    plt.ylabel('t-SNE Component 2')

    from matplotlib.lines import Line2D

    handles = [
        Line2D([], [], marker='o', color=colors[i], linestyle='',
            markersize=10, label=uniques[i])
        for i in range(n_classes)
    ]

    plt.legend(handles=handles, title="Categories",
            bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

# feature selection: Lasso (L1 regularization)
# 12k features is a lot...
# the code below is wrong because y_train has 65 columns rather than 1
# def lasso_feature_selection(X_train, X_test, feature_names, y_train, alpha=0.01):
#     """
#     Only keep features with non-zero Lasso coefficients
#     """
#     from sklearn.linear_model import Lasso
#     import numpy as np
    
#     print(f"Original feature count: {X_train.shape[1]}")
    
#     # Fit Lasso
#     lasso = Lasso(alpha=alpha, random_state=42, max_iter=10000)
#     lasso.fit(X_train, y_train)
    
#     # Get non-zero coefficients
#     non_zero_mask = lasso.coef_ != 0
#     n_selected = np.sum(non_zero_mask)
    
#     print(f"Features with non-zero coefficients: {n_selected}")
    
#     # Apply selection manually
#     X_train_selected = X_train[:, non_zero_mask]
#     X_test_selected = X_test[:, non_zero_mask]
    
#     selected_features = [feature_names[i] for i in range(len(feature_names)) if non_zero_mask[i]]
    
#     print(f"Selected features: {X_train_selected.shape[1]} ({n_selected/X_train.shape[1]:.1%} kept)")
    
#     return {
#         'X_train_selected': X_train_selected,
#         'X_test_selected': X_test_selected,
#         'selected_features': selected_features,
#         'non_zero_mask': non_zero_mask
#     }

# new approach: run Lasso per label and combine the masks, keeping features that are useful for at least one label 
from sklearn.linear_model import Lasso

def lasso_feature_selection(X_train, X_test, feature_names, y_train, alpha=0.01):
    n_labels = y_train.shape[1]
    n_features = X_train.shape[1]
    combined_mask = np.zeros(n_features, dtype=bool)
    for label_idx in range(n_labels):
        lasso = Lasso(alpha=alpha, random_state=42, max_iter=10000)
        lasso.fit(X_train, y_train.iloc[:, label_idx])
        combined_mask |= (lasso.coef_ != 0)
    X_train_selected = X_train.loc[:, combined_mask]
    X_test_selected  = X_test.loc[:, combined_mask]
    selected_features = [feature_names[i] for i in range(len(feature_names)) if combined_mask[i]]
    print(f"Selected {combined_mask.sum()} features out of {n_features}")
    return {
        "X_train_selected": X_train_selected,
        "X_test_selected": X_test_selected,
        "selected_features": selected_features,
        "mask": combined_mask
    }


# Select features using Lasso (we aren't using Lasso to classify yet)
feature_names = X_df.columns.tolist()  # Get feature names

X_df_train = X_df.iloc[train_idx]
X_df_test = X_df.iloc[test_idx]

selection_results = lasso_feature_selection(
    X_train=X_df_train,
    X_test=X_df_test,
    feature_names=feature_names,
    y_train=y_train,
    alpha=0.01
)
# Selected 13 features out of 12781
# Lasso is clearly not the right tool for this kind of data
# will try other methods

# def feature_selection():
# Strategy: Use chi2 to filter down first, then use RF for refinement
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.preprocessing import MultiLabelBinarizer

# If your y is already in multi-label format (0/1 matrix)
# If not, binarize it:
# mlb = MultiLabelBinarizer()
# y_train_bin = mlb.fit_transform(y_train)
# y_test_bin = mlb.transform(y_test)

# For chi2, you need to compute scores per label
# Approach: Aggregate chi2 scores across all labels
from sklearn.feature_selection import chi2
import numpy as np

# def select_features_chi2_multilabel(X, y, k=1000):
#     """
#     Select top k features using chi2 for multi-label data.
#     Uses maximum chi2 score across all labels for each feature.
#     """
#     n_features = X.shape[1]
#     chi2_scores = np.zeros(n_features)
#     # For each label, compute chi2 scores
#     for label_idx in range(y.shape[1]):
#         chi2_vals, p_vals = chi2(X, y[:, label_idx])
#         chi2_scores = np.maximum(chi2_scores, chi2_vals)  # Take max across labels
#     # Get top k features
#     top_indices = np.argsort(chi2_scores)[-k:][::-1]
#     return top_indices, chi2_scores

# # Select top 1000 features (adjust k based on your needs)
# k_features = 1000
# selected_indices_chi2, chi2_scores = select_features_chi2_multilabel(X_train, y_train, k=k_features)
# # Transform your data
# X_train_chi2 = X_train[:, selected_indices_chi2]
# X_test_chi2 = X_test[:, selected_indices_chi2]
# print(f"Original features: {X_train.shape[1]}")
# print(f"Selected features: {X_train_chi2.shape[1]}")

from sklearn.feature_selection import SelectKBest, chi2
from sklearn.multioutput import MultiOutputClassifier

# If you want to use the average chi2 score across labels
def multi_label_chi2(X, y):
    """Compute chi2 for multi-label by averaging scores."""
    scores = []
    for i in range(y.shape[1]):
        chi2_vals, _ = chi2(X, y[:, i])
        scores.append(chi2_vals)
    return np.mean(scores, axis=0), None

# Using SelectKBest with custom scoring
selector_chi2 = SelectKBest(score_func=multi_label_chi2, k=1000)
X_train_chi2 = selector_chi2.fit_transform(X_train, y_train)
X_test_chi2 = selector_chi2.transform(X_test)

# Get selected feature indices
selected_indices_chi2 = selector_chi2.get_support(indices=True)

from sklearn.ensemble import RandomForestClassifier
import numpy as np

def select_features_rf_multilabel(X, y, k=1000, n_estimators=100):
    """
    Select top k features using Random Forest for multi-label data.
    """
    # For multi-label, use RandomForest with multi-output
    # This trains one tree per label
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=42,
        n_jobs=-1,
        class_weight='balanced'  # Good for imbalanced labels
    )
    # Fit the model
    rf.fit(X, y)
    # Get feature importances (average across all trees)
    # For multi-label, we can average importance across labels
    feature_importances = rf.feature_importances_
    # Get top k features
    top_indices = np.argsort(feature_importances)[-k:][::-1]
    return top_indices, feature_importances, rf

# Select top 1000 features
k_features = 1000
selected_indices_rf, rf_importances, rf_model = select_features_rf_multilabel(
    X_train, y_train, k=k_features
)

# Transform your data
X_train_rf = X_train[:, selected_indices_rf]
X_test_rf = X_test[:, selected_indices_rf]

print(f"Original features: {X_train.shape[1]}")
print(f"Selected features: {X_train_rf.shape[1]}")

def hybrid_feature_selection(X, y, chi2_k=3000, rf_k=1000):
    """
    Two-stage feature selection:
    1. Use chi2 to reduce to chi2_k features
    2. Use Random Forest to select top rf_k from those
    """
    # Stage 1: Chi2 selection
    selector_chi2 = SelectKBest(score_func=multi_label_chi2, k=chi2_k)
    X_chi2 = selector_chi2.fit_transform(X, y)
    chi2_indices = selector_chi2.get_support(indices=True)
    # Stage 2: Random Forest on reduced set
    rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
    rf.fit(X_chi2, y)
    # Get top rf_k features from the reduced set
    rf_importances = rf.feature_importances_
    rf_top_local = np.argsort(rf_importances)[-rf_k:][::-1]
    # Map back to original indices
    final_indices = chi2_indices[rf_top_local]
    return final_indices

# Apply hybrid selection
selected_indices_hybrid = hybrid_feature_selection(
    X_train, y_train, chi2_k=3000, rf_k=1000
)
X_train_hybrid = X_train[:, selected_indices_hybrid]
X_test_hybrid = X_test[:, selected_indices_hybrid]

from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import f1_score, accuracy_score

# def evaluate_features(X_train, X_test, y_train, y_test, method_name):
#     """Quick evaluation of feature set."""
#     clf = MultiOutputClassifier(
#         LogisticRegression(max_iter=1000, random_state=42),
#         n_jobs=-1
#     )
#     clf.fit(X_train, y_train)
#     y_pred = clf.predict(X_test)
#     # Use macro F1 for multi-label
#     f1 = f1_score(y_test, y_pred, average='macro')
#     print(f"{method_name} - Macro F1: {f1:.4f}")

#     return f1

from sklearn.metrics import f1_score, jaccard_score
import numpy as np

def evaluate_features(X_train, X_test, y_train, y_test, method_name):
    """Quick evaluation of feature set."""
    clf = MultiOutputClassifier(
        LogisticRegression(max_iter=1000, random_state=42),
        n_jobs=-1
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    # Compute macro F1
    f1 = f1_score(y_test, y_pred, average='macro')
    # Compute per-sample Jaccard by treating each sample as a multi-label problem
    # jaccard_score with average='samples' computes the average Jaccard across samples
    avg_jaccard = jaccard_score(y_test, y_pred, average='samples', zero_division=0)
    print(f"{method_name} - Macro F1: {f1:.4f}, Avg Jaccard: {avg_jaccard:.4f}")
    return f1, avg_jaccard

# Compare different selections
results = {}
results['original'] = evaluate_features(X_train, X_test, y_train, y_test, 'Original')
results['chi2'] = evaluate_features(X_train_chi2, X_test_chi2, y_train, y_test, 'Chi2')
results['rf'] = evaluate_features(X_train_rf, X_test_rf, y_train, y_test, 'Random Forest')
results['hybrid'] = evaluate_features(X_train_hybrid, X_test_hybrid, y_train, y_test, 'Hybrid')

print("Scores:")
for k, v in results.items():
    print(f"{k}: Macro F1 {v[0]:.4f}, Average Jaccard {v[1]:.4f}")
# original: Macro F1 0.5504, Average Jaccard 0.8141
# chi2: Macro F1 0.5415, Average Jaccard 0.7855
# rf: Macro F1 0.4717, Average Jaccard 0.7249
# hybrid: Macro F1 0.5350, Average Jaccard 0.8000
# based on the Macro F1 scores, looks like chi2 is the clear winner, but a macro F1 score of 0.5 out of 1 is borderline.
# probably struggles to predict the rare classes. Hopefully models other than logistic regression do better.
# Then again, we don't necessarily expect the model to predict perfectly on every multi-label classification problem since there are multiple labels to assign or not
# So maybe the Jaccard is a better index to use here.
# Hm, the hybrid approach has a lower F1 but higher Jaccard as compared to the chi2. Actually, the hybrid did really well if we base it on the Jaccard.
# Considering there were 12k features to start with, 
# For feature selection we're looking for a method that has F1 and/or Jaccard close to (but inevitably slightly lower than) the original
# Btw the reason the Jaccard for the original slightly differs from the one printed in the first Jaccard score test above is that I used n_to_test to limit the number of samples to evaluate.

# ok, in the real analysis you can apply the hybrid feature selection function and use X_train_hybrid and X_test_hybrid, rather than their supersets X_train and X_test, to train and test the models.

# k_values = [500, 1000, 2000, 3000, 5000]
k_values = [1000, 1250, 1500, 1750, 2000]
for i in range(len(k_values)):
    chi2_k = k_values[i]
    rf_k = k_values[i] # I used to use a different value for this, but for a fair comparison let's just use the same as chi2_k
    print(chi2_k, rf_k)
    #
    selector_chi2 = SelectKBest(score_func=multi_label_chi2, k=chi2_k)
    X_train_chi2 = selector_chi2.fit_transform(X_train, y_train)
    X_test_chi2 = selector_chi2.transform(X_test)
    #
    selected_indices_rf, rf_importances, rf_model = select_features_rf_multilabel(
        X_train, y_train, k=rf_k
    )
    X_train_rf = X_train[:, selected_indices_rf]
    X_test_rf = X_test[:, selected_indices_rf]
    #
    selected_indices_hybrid = hybrid_feature_selection(
        X_train, y_train, chi2_k=2*chi2_k, rf_k=rf_k
    )
    X_train_hybrid = X_train[:, selected_indices_hybrid]
    X_test_hybrid = X_test[:, selected_indices_hybrid]
    #
    results = {}
    results['original'] = evaluate_features(X_train, X_test, y_train, y_test, 'Original')
    results['chi2'] = evaluate_features(X_train_chi2, X_test_chi2, y_train, y_test, 'Chi2')
    results['rf'] = evaluate_features(X_train_rf, X_test_rf, y_train, y_test, 'Random Forest')
    results['hybrid'] = evaluate_features(X_train_hybrid, X_test_hybrid, y_train, y_test, 'Hybrid')
    #
    # print("Scores:")
    # for k, v in results.items():
    #     print(f"{k}: Macro F1 {v[0]:.4f}, Average Jaccard {v[1]:.4f}")

# 500 500
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5274, Avg Jaccard: 0.7568
# Random Forest - Macro F1: 0.3880, Avg Jaccard: 0.6446
# Hybrid - Macro F1: 0.5131, Avg Jaccard: 0.7806
# 1000 1000
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5415, Avg Jaccard: 0.7855
# Random Forest - Macro F1: 0.4717, Avg Jaccard: 0.7249
# Hybrid - Macro F1: 0.5308, Avg Jaccard: 0.7970
# 2000 2000
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5462, Avg Jaccard: 0.7992
# Random Forest - Macro F1: 0.5276, Avg Jaccard: 0.7774
# Hybrid - Macro F1: 0.5488, Avg Jaccard: 0.8083
# 3000 3000
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5490, Avg Jaccard: 0.8056
# Random Forest - Macro F1: 0.5289, Avg Jaccard: 0.7837
# Hybrid - Macro F1: 0.5506, Avg Jaccard: 0.8107
# 5000 5000
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5500, Avg Jaccard: 0.8114
# Random Forest - Macro F1: 0.5384, Avg Jaccard: 0.7926
# Hybrid - Macro F1: 0.5495, Avg Jaccard: 0.8122
# based on this: hybrid is definitely the best, but what k should I use? Check for diminishing returns.

# honing in on the k=1000 to k=2000 range: eh, I still feel like k=2000 is appropriate.
# 1000 1000
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5415, Avg Jaccard: 0.7855
# Random Forest - Macro F1: 0.4717, Avg Jaccard: 0.7249
# Hybrid - Macro F1: 0.5308, Avg Jaccard: 0.7970
# 1250 1250
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5467, Avg Jaccard: 0.7909
# Random Forest - Macro F1: 0.4930, Avg Jaccard: 0.7556
# Hybrid - Macro F1: 0.5479, Avg Jaccard: 0.8009
# 1500 1500
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5450, Avg Jaccard: 0.7929
# Random Forest - Macro F1: 0.5115, Avg Jaccard: 0.7629
# Hybrid - Macro F1: 0.5498, Avg Jaccard: 0.8032
# 1750 1750
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5464, Avg Jaccard: 0.7969
# Random Forest - Macro F1: 0.5223, Avg Jaccard: 0.7705
# Hybrid - Macro F1: 0.5504, Avg Jaccard: 0.8058
# 2000 2000
# Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
# Chi2 - Macro F1: 0.5462, Avg Jaccard: 0.7992
# Random Forest - Macro F1: 0.5276, Avg Jaccard: 0.7774
# Hybrid - Macro F1: 0.5488, Avg Jaccard: 0.8083

import matplotlib.pyplot as plt

def plot_feature_importance(importances, title, top_n=20):
    """Plot top N feature importances."""
    top_indices = np.argsort(importances)[-top_n:]
    top_importances = importances[top_indices]
    plt.figure(figsize=(10, 6))
    plt.barh(range(top_n), top_importances)
    plt.yticks(range(top_n), [f'Feature {i}' for i in top_indices])
    plt.xlabel('Importance Score')
    plt.title(title)
    plt.tight_layout()
    plt.show()

# Plot RF feature importances
plot_feature_importance(rf_importances, 'Random Forest Feature Importance (Top 20)')


# from sklearn.utils import resample

# def stability_analysis(X, y, n_iterations=10, k=1000):
#     """Check stability of feature selection across bootstrap samples."""
#     selected_sets = []
#     for i in range(n_iterations):
#         # Bootstrap sample
#         X_boot, y_boot = resample(X, y, random_state=i)
#         indices, _ = select_features_rf_multilabel(X_boot, y_boot, k=k)
#         selected_sets.append(set(indices))
#     # Calculate Jaccard similarity between all pairs
#     similarities = []
#     for i in range(len(selected_sets)):
#         for j in range(i+1, len(selected_sets)):
#             jaccard = len(selected_sets[i] & selected_sets[j]) / len(selected_sets[i] | selected_sets[j])
#             similarities.append(jaccard)
#     print(f"Average Jaccard similarity: {np.mean(similarities):.3f}")
#     print(f"Std: {np.std(similarities):.3f}")

# stability_analysis(X_train, y_train)

from sklearn.metrics import f1_score, jaccard_score, accuracy_score, hamming_loss
import numpy as np
import pandas as pd

def evaluate_features_comprehensive(X_train, X_test, y_train, y_test, method_name):
    """Comprehensive evaluation with multiple metrics."""
    clf = MultiOutputClassifier(
        LogisticRegression(max_iter=1000, random_state=42),
        n_jobs=-1
    )
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    
    # Multiple metrics
    metrics = {}
    
    # 1. Jaccard (samples average) - YOUR PREFERRED METRIC
    metrics['jaccard_samples'] = jaccard_score(y_test, y_pred, average='samples', zero_division=0)
    
    # 2. Macro F1 (average across labels)
    metrics['f1_macro'] = f1_score(y_test, y_pred, average='macro', zero_division=0)
    
    # 3. Micro F1 (global counts)
    metrics['f1_micro'] = f1_score(y_test, y_pred, average='micro', zero_division=0)
    
    # 4. Per-label F1 (to see rare label performance)
    per_label_f1 = f1_score(y_test, y_pred, average=None, zero_division=0)
    
    # 5. Hamming Loss (fraction of incorrect labels)
    metrics['hamming_loss'] = hamming_loss(y_test, y_pred)
    
    # 6. Exact Match Ratio (strict accuracy)
    metrics['exact_match'] = accuracy_score(y_test, y_pred)
    
    # Print results
    print(f"\n{method_name} Results:")
    print("-" * 40)
    print(f"Jaccard (samples avg): {metrics['jaccard_samples']:.4f}")
    print(f"F1 Macro:              {metrics['f1_macro']:.4f}")
    print(f"F1 Micro:              {metrics['f1_micro']:.4f}")
    print(f"Hamming Loss:          {metrics['hamming_loss']:.4f}")
    print(f"Exact Match:           {metrics['exact_match']:.4f}")
    
    # Show per-label performance
    print("\nPer-label F1 scores:")
    for i, f1 in enumerate(per_label_f1):
        label_positive_count = y_test[:, i].sum()
        print(f"  Label {i}: F1={f1:.4f} (positives: {label_positive_count})")
    
    return metrics