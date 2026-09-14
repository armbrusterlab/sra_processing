# may clean up and encapsulate as functions later, but for now, CD to a new directory to avoid overwriting old results 

import pandas as pd
import numpy as np
import joblib

#df = pd.read_csv(r"C:\Users\achro\OneDrive\Desktop\CMU\Spring 2025\Armbruster Lab research\bacdiveReformat_2026-08-12.tsv",
                 #sep="\t") # processed in wrangle_bacdive.py to add joined_1_2 column
df = pd.read_csv("/home/kcw2/data/testing/bacdive_model/bacdiveReformat_2026-08-12.tsv",
                 sep="\t") # processed in wrangle_bacdive.py to add joined_1_2 column
                 
print(f"Original length of df: {len(df)}") # 63296

# row-wise filtering
df = df[df["Isolation source"].isna() == False]
print(f"Length of df after removing rows with no isolation source: {len(df)}") # 58625

df = df[df["joined_1_2"].isna() == False]
print(f"Length of df after removing rows with no category tags: {len(df)}") # 43314

# reconvert joined_1_2 from string to list
df["joined_1_2"] = [
    row["joined_1_2"].split("###")
    for _, row in df.iterrows()
]

terms = {} # frequency dictionary
for t in df["joined_1_2"]:
    for term in t:
        terms[term] = terms.get(term, 0) + 1

# In bacdive_eda.py, I found that there were 117 terms. 30 terms appeared only once.
# For training, I will only keep terms that appear at least 20 times.
minSize = 20
terms_blacklist = set([k for k in terms.keys() if terms[k] < minSize])
terms_blacklist = terms_blacklist.union(set([k for k in terms.keys() if "Host@@@" in k])) # decided to filter out "Host" category due to its ambiguity
terms_blacklist.add("Infection@@@Patient") # this term in particular seems to be ambiguous and difficult to categorize
terms_keep = sorted(list(set(terms.keys()) - terms_blacklist)) # Fixing a MAJOR bug- if this is a set, it's impossible to actually reconstitute the predicted terms from the prediction columns

# will need terms_keep later in order to decipher the model's y_pred output
# For now, create the output dir in the current directory
import os
output_folder = "models"
os.makedirs(output_folder, exist_ok=True)

joblib.dump(terms_keep, 'models/y_colnames.joblib')

# Filter out infrequent terms from the joined_1_2 column
df["joined_1_2"] = [
    list(set(row["joined_1_2"]) - terms_blacklist)
    for _, row in df.iterrows()
]

# Some rows may now have empty lists after removing terms that are too rare, so remove these rows
df = df[df["joined_1_2"].str.len() > 0] # even though the items are lists, pandas's str.len() function can get list lengths
print(f"Length of df after removing extremely rare terms: {len(df)}") # 43283; compared to 43314, didn't remove too much

# Convert the response column, joined_1_2, into binary yes/no columns, one per label
y = pd.DataFrame()

for label in terms_keep:
    # print(f"Processing {label}...")
    y[label] = [1 if label in t else 0 for t in df["joined_1_2"]]

# To make the predictor columns (the features), vectorize the isolation source text data
# https://www.geeksforgeeks.org/nlp/text-classification-using-scikit-learn-in-nlp/
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction import text

# add bacteria names as stopwords
# Bacteria name list obtained from https://app.gideononline.com/az/pathogens/bacteria
with open("/home/kcw2/metadata-magnet/scripts/modeling/bacteria_stopwords.txt", "r", errors="ignore") as file:
   bacteria_str = file.read()
bacteria = set([s.lower() for s in bacteria_str.split() if not s.isnumeric()]) # separate genus and species names; split by any whitespace
bacteria.remove("(commensal)")

# add country names as stopwords
# country list obtained from https://gist.github.com/dariusz-wozniak/656f2f9070b4205c5009716f05c94067
with open("/home/kcw2/metadata-magnet/scripts/modeling/country_stopwords.txt", "r", errors="ignore") as file:
   country_str = file.read()
countries = set([s.lower() for s in country_str.split() if not s.isnumeric()]) # separate genus and species names; split by any whitespace
countries.remove("&")
countries.remove("island")
countries.remove("islands")
countries.remove("ocean")
countries.remove("city")
countries.remove("sandwich")

# add other stopwords
countries.add("collected") # the phrase "not collected" was being marked with a variety of categories for some reason
countries.update(["missing", "unknown", "metagenome", "not", "provided"]) # these words don't seem to cause issues, but they commonly occur and are uninformative, so I'll add them in case
countries.update("metagenomes available from the Sequence Read Archive".lower().split()) # this has the potential to be problematic too
countries.add("sterile") # the models associate this with laboratory, but it might refer to a sterile body site, so it's not informative

# add integers as stopwords
integer_stopwords = set([str(i) for i in range(100)]) # manual inspection of selected features suggests that 2 digits are enough

custom_stopwords = list(text.ENGLISH_STOP_WORDS.union(bacteria).union(countries).union(integer_stopwords))

# Initialize TF-IDF Vectorizer
vectorizer = TfidfVectorizer(stop_words=custom_stopwords, max_df=0.7) # first-pass "feature selection"; remove words that appear in >70% of the records since they're not helpful
# I could use min_df but opted not to, as many groups are very small and this might remove legitimate information

# vectorizer = TfidfVectorizer(stop_words='english', max_df=0.7)

# Transform the text data to feature vectors
# first, convert the ### delimiter in the isolation source strings into spaces so that the words on either side of the delimiter will be processed properly
# also remove the <I> and </I> for italicization
iso_source = [
    row["Isolation source"].replace("###", " ")
        .replace("<I>", "").replace("</I>", "")
        .replace("<i>", "").replace("</i>", "")
        .replace("-", " ")
        .replace("(", "").replace(")", "")
        .replace('"', '')
        .lower()
    for _, row in df.iterrows()
]
X = vectorizer.fit_transform(iso_source)
print(f"Dimensions of X: {X.shape}") # (43283, 12781); there are over 12k features to begin with.

# Save the fitted vectorizer for later use
joblib.dump(vectorizer, 'models/tfidf_vectorizer.joblib')

# I want to know the indices corresponding to the train-test split, so I will split the indices rather than splitting the df directly
from sklearn.model_selection import train_test_split

indices = np.arange(len(df))
train_idx, test_idx = train_test_split(
    indices, test_size=0.2, random_state=42
)

X_train = X[train_idx]
X_test  = X[test_idx]
y_train = y.iloc[train_idx]
y_test  = y.iloc[test_idx]

# Feature selection: during EDA I determined that it works best with a hybrid of chi square and random forest approaches
# and that by narrowing it down to the top 2000 features we can retain much of the signal while removing the noise
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.ensemble import RandomForestClassifier

# use the average chi2 score across labels
def multi_label_chi2(X, y):
    """Compute chi2 for multi-label by averaging scores."""
    scores = []
    for i in range(y.shape[1]):
        chi2_vals, _ = chi2(X, y[:, i])
        scores.append(chi2_vals)
    return np.mean(scores, axis=0), None

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

from sklearn.metrics import f1_score, jaccard_score
from sklearn.linear_model import LogisticRegression # below, it's used as a quick way to evaluate model performance
from sklearn.multioutput import MultiOutputClassifier

def evaluate_features(X_train, X_test, y_train, y_test, method_name):
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

# One last bit of EDA: I decided on an ending k (rf_k) of 2000, but for the starting k (chi2_k), should it be larger or smaller?
def eda():
    # Apply hybrid selection
    selected_indices_hybrid = hybrid_feature_selection(
        X_train, y_train, chi2_k=4000, rf_k=2000
    )
    X_train_hybrid = X_train[:, selected_indices_hybrid]
    X_test_hybrid = X_test[:, selected_indices_hybrid]

    selected_indices_hybrid = hybrid_feature_selection(
        X_train, y_train, chi2_k=6000, rf_k=2000
    )
    X_train_hybrid2 = X_train[:, selected_indices_hybrid]
    X_test_hybrid2 = X_test[:, selected_indices_hybrid]

    selected_indices_hybrid = hybrid_feature_selection(
        X_train, y_train, chi2_k=8000, rf_k=2000
    )
    X_train_hybrid3 = X_train[:, selected_indices_hybrid]
    X_test_hybrid3 = X_test[:, selected_indices_hybrid]

    # Compare selections
    results = {}
    results['original'] = evaluate_features(X_train, X_test, y_train, y_test, 'Original') # Original - Macro F1: 0.5504, Avg Jaccard: 0.8141
    results['hybrid4000'] = evaluate_features(X_train_hybrid, X_test_hybrid, y_train, y_test, 'Hybrid: chi2_k=4000') # Hybrid: chi2_k=4000 - Macro F1: 0.5488, Avg Jaccard: 0.8083
    results['hybrid6000'] = evaluate_features(X_train_hybrid2, X_test_hybrid2, y_train, y_test, 'Hybrid: chi2_k=6000') # Hybrid: chi2_k=6000 - Macro F1: 0.5492, Avg Jaccard: 0.8069
    results['hybrid8000'] = evaluate_features(X_train_hybrid3, X_test_hybrid3, y_train, y_test, 'Hybrid: chi2_k=8000') # Hybrid: chi2_k=8000 - Macro F1: 0.5509, Avg Jaccard: 0.8090

# Since the difference between the trials is small (using Jaccard index as the main scoring criteria), I will just go with chi2_k=4000
chi2_k = 4000 # X_train.shape[1]
rf_k = 2000 # X_train.shape[1] #1000
selected_indices_hybrid = hybrid_feature_selection(
   X_train, y_train, chi2_k=chi2_k, rf_k=rf_k
)
X_train_hybrid = X_train[:, selected_indices_hybrid]
X_test_hybrid = X_test[:, selected_indices_hybrid]
# These dataframes contain the top 2000 features (or more generally speaking, rf_k number of features).
# >>> X_train_hybrid.shape  
# (34626, 2000)
# >>> X_test_hybrid.shape 
# (8657, 2000)

# # minimally disruptive (in terms of pipeline) way to skip feature selection: just select all indices...
# selected_indices_hybrid = hybrid_feature_selection(
#    X_train, y_train, chi2_k=X_train.shape[1], rf_k=X_train.shape[1]
# )
# X_train_hybrid = X_train[:, selected_indices_hybrid]
# X_test_hybrid = X_test[:, selected_indices_hybrid]

X_train_hybrid.shape 
X_test_hybrid.shape 

# Also need this later when predicting upon new data
joblib.dump(selected_indices_hybrid, 'models/feature_selection_indices.joblib')

# No scaling necessary because all of the features are binary.

# Train and assess the models:
# The code below was generated with DeepSeek.
import optuna
import time
import logging
import pickle
import os
from datetime import datetime
from sklearn.model_selection import cross_val_score, StratifiedKFold, KFold
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.linear_model import LogisticRegression, Lasso
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import f1_score, jaccard_score, accuracy_score, precision_score
from sklearn.preprocessing import LabelEncoder
from sklearn.base import clone
import warnings
warnings.filterwarnings('ignore')

def build_and_tune_models(X_train, y_train, X_test, y_test, n_trials=10, save_dir='models'):
    """
    Build and tune multiple models using Optuna with multi-label support.
    
    Parameters:
    -----------
    X_train, y_train : Training data
    X_test, y_test : Test data
    n_trials : Number of Optuna trials per model
    save_dir : Directory to save the best model
    
    Returns:
    --------
    dict : Contains all trained models and results
    """
    
    # Create save directory if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('bacdive_tuning.log'),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    
    logger.info("="*80)
    logger.info(f"Starting model tuning at {datetime.now()}")
    logger.info(f"Training data shape: {X_train.shape}")
    logger.info(f"Number of labels: {y_train.shape[1]}")
    logger.info(f"Best model will be saved to: {save_dir}")
    logger.info("="*80)
    
    results = {
        'models': {},
        'best_params': {},
        'best_scores': {},
        'tuning_times': {},
        'cv_scores': {}
    }
    
    # Custom scorer for cross-validation (Jaccard)
    def jaccard_scorer(estimator, X, y):
        """Custom scorer for multi-label Jaccard score."""
        y_pred = estimator.predict(X)
        return jaccard_score(y, y_pred, average='samples', zero_division=0)
    
    # Define parameter distributions for each model
    param_distributions = {
        'random_forest': {
            'estimator__n_estimators': optuna.distributions.IntDistribution(50, 300),
            'estimator__max_depth': optuna.distributions.CategoricalDistribution([5, 10, 15, 20, None]),
            'estimator__min_samples_split': optuna.distributions.IntDistribution(2, 20),
            'estimator__min_samples_leaf': optuna.distributions.IntDistribution(1, 10),
            'estimator__max_features': optuna.distributions.CategoricalDistribution(['sqrt', 'log2', None]),
            'estimator__random_state': optuna.distributions.CategoricalDistribution([42])
        },
        'gradient_boosting': {
            'estimator__n_estimators': optuna.distributions.IntDistribution(50, 200),
            'estimator__learning_rate': optuna.distributions.FloatDistribution(0.01, 0.3),
            'estimator__max_depth': optuna.distributions.IntDistribution(3, 10),
            'estimator__min_samples_leaf': optuna.distributions.IntDistribution(10, 50),
            'estimator__random_state': optuna.distributions.CategoricalDistribution([42])
        },
        'lasso': {
            'estimator__C': optuna.distributions.FloatDistribution(0.001, 1.0),
            'estimator__penalty': optuna.distributions.CategoricalDistribution(['l1']),
            'estimator__solver': optuna.distributions.CategoricalDistribution(['liblinear', 'saga']),
            'estimator__max_iter': optuna.distributions.CategoricalDistribution([1000]),
            'estimator__random_state': optuna.distributions.CategoricalDistribution([42])
        },
        'logistic_regression': {
            'estimator__penalty': optuna.distributions.CategoricalDistribution(['l1', 'l2', 'elasticnet', None]),
            'estimator__C': optuna.distributions.FloatDistribution(0.0001, 10000, log=True),
            'estimator__solver': optuna.distributions.CategoricalDistribution(['lbfgs', 'liblinear', 'saga']),
            'estimator__max_iter': optuna.distributions.CategoricalDistribution([100, 1000, 2500]),
            'estimator__random_state': optuna.distributions.CategoricalDistribution([42])
        }
    }
    
    # Define model creation functions
    def create_random_forest(params):
        rf = RandomForestClassifier(
            n_estimators=params['estimator__n_estimators'],
            max_depth=params['estimator__max_depth'],
            min_samples_split=params['estimator__min_samples_split'],
            min_samples_leaf=params['estimator__min_samples_leaf'],
            max_features=params['estimator__max_features'],
            random_state=42,
            n_jobs=-1
        )
        return MultiOutputClassifier(rf, n_jobs=-1)
    
    def create_gradient_boosting(params):
        gb = GradientBoostingClassifier(
            n_estimators=params['estimator__n_estimators'],
            learning_rate=params['estimator__learning_rate'],
            max_depth=params['estimator__max_depth'],
            min_samples_leaf=params['estimator__min_samples_leaf'],
            random_state=42
        )
        return MultiOutputClassifier(gb, n_jobs=-1)
    
    def create_lasso(params):
        # Using LogisticRegression with l1 penalty as an alternative to Lasso
        lr = LogisticRegression(
            C=params['estimator__C'],
            penalty='l1',
            solver=params['estimator__solver'],
            max_iter=params['estimator__max_iter'],
            random_state=42
        )
        return MultiOutputClassifier(lr, n_jobs=-1)
    
    def create_logistic_regression(params):
        """Create logistic regression model with parameter validation."""
        
        # Get the penalty
        penalty = params.get('estimator__penalty')
        solver = params.get('estimator__solver')
        
        # Validate and fix parameter combinations
        if penalty == 'elasticnet':
            # elasticnet only works with saga
            solver = 'saga'
            # Must have l1_ratio for elasticnet
            if 'estimator__l1_ratio' not in params:
                params['estimator__l1_ratio'] = 0.5
        elif penalty == 'l1':
            # l1 works with liblinear or saga
            if solver not in ['liblinear', 'saga']:
                solver = 'liblinear'  # default to liblinear
        elif penalty == 'l2':
            # l2 works with most solvers
            if solver not in ['lbfgs', 'liblinear', 'saga']:
                solver = 'lbfgs'  # default to lbfgs
        else:  # penalty == None
            # None penalty works with lbfgs or saga (not liblinear)
            if solver not in ['lbfgs', 'saga']:
                solver = 'lbfgs'  # default to lbfgs
        
        # Build the model parameters
        lr_params = {
            'penalty': penalty,
            'C': params.get('estimator__C', 1.0),
            'solver': solver,
            'max_iter': params.get('estimator__max_iter', 1000),
            'random_state': 42
        }
        
        # Add l1_ratio if using elasticnet
        if penalty == 'elasticnet':
            lr_params['l1_ratio'] = params.get('estimator__l1_ratio', 0.5)
        
        return MultiOutputClassifier(LogisticRegression(**lr_params), n_jobs=-1)
    
    def objective(trial, model_name, create_model_func, param_dist):
        """Objective function for Optuna optimization."""
        # Sample parameters
        params = {}
        for param_name, dist in param_dist.items():
            if isinstance(dist, optuna.distributions.CategoricalDistribution):
                params[param_name] = trial.suggest_categorical(param_name, dist.choices)
            elif isinstance(dist, optuna.distributions.IntDistribution):
                params[param_name] = trial.suggest_int(param_name, dist.low, dist.high)
            elif isinstance(dist, optuna.distributions.FloatDistribution):
                params[param_name] = trial.suggest_float(param_name, dist.low, dist.high, log=dist.log)
        
        # Create model
        try:
            print(params)
            model = create_model_func(params)
            
            # Cross-validation with Jaccard score
            cv = KFold(n_splits=5, shuffle=True, random_state=42)
            
            # Compute cross-validation scores
            scores = cross_val_score(
                model, X_train, y_train,
                cv=cv,
                # scoring=jaccard_scorer,
                # scoring='jaccard_weighted',
                # scoring='jaccard_macro',
                # scoring='f1_macro',
                # scoring='f1_weighted',
                scoring='precision_macro',
                # scoring='precision_weighted',
                n_jobs=-1,
                error_score='raise'
            )
            
            return scores.mean()
            
        except Exception as e:
            logger.warning(f"Trial failed: {str(e)}")
            return -1.0
    
    # Optimize each model type
    # model_configs = [
    #     ('random_forest', create_random_forest, param_distributions['random_forest']),
    #     ('gradient_boosting', create_gradient_boosting, param_distributions['gradient_boosting']),
    #     ('lasso', create_lasso, param_distributions['lasso']),
    #     ('logistic_regression', create_logistic_regression, param_distributions['logistic_regression'])
    # ]

    # version that only runs trials for logistic regression
    model_configs = [
        ('logistic_regression', create_logistic_regression, param_distributions['logistic_regression'])
    ]
    
    # Store trained base models for stacking
    base_models = []
    
    for model_name, create_func, param_dist in model_configs:
        logger.info(f"\n{'='*60}")
        logger.info(f"Optimizing {model_name.upper()}")
        logger.info(f"{'='*60}")
        
        start_time = time.time()
        
        # Create Optuna study
        study = optuna.create_study(
            direction='maximize',
            sampler=optuna.samplers.RandomSampler(seed=42)
        )
        
        # Optimize
        study.optimize(
            lambda trial: objective(trial, model_name, create_func, param_dist),
            n_trials=n_trials,
            show_progress_bar=True
        )
        
        tuning_time = time.time() - start_time
        
        # Get best parameters and score
        best_params = study.best_params
        best_score = study.best_value
        
        # Train best model on full training data
        logger.info(f"\nTraining best {model_name} on full training data...")
        
        # Convert params back to original format for create_func
        create_params = {k: v for k, v in best_params.items()}
        best_model = create_func(create_params)
        best_model.fit(X_train, y_train)
        
        # Evaluate on test set
        y_pred = best_model.predict(X_test)
        f1 = f1_score(y_test, y_pred, average='macro', zero_division=0)
        jaccard = jaccard_score(y_test, y_pred, average='samples', zero_division=0)
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred, average='macro', zero_division=0)
        
        scores = {
            'cv_score': best_score,
            'test_f1': f1,
            'test_jaccard': jaccard,
            'test_accuracy': accuracy,
            'test_precision': precision
        }
        
        # Store results (but not the model itself to save memory)
        results['models'][model_name] = best_model
        results['best_params'][model_name] = best_params
        results['best_scores'][model_name] = scores
        results['tuning_times'][model_name] = tuning_time
        results['cv_scores'][model_name] = study.trials_dataframe()['value'].values
        
        # Store best base model for stacking (if not lasso)
        if model_name != 'lasso':
            if isinstance(best_model, MultiOutputClassifier):
                base_models.append((model_name, best_model.estimator))
            else:
                base_models.append((model_name, best_model))
        
        # Print results immediately
        logger.info(f"\n{model_name.upper()} Results:")
        logger.info(f"  Best CV score: {best_score:.4f}")
        logger.info(f"  Test F1 (macro): {f1:.4f}")
        logger.info(f"  Test Jaccard: {jaccard:.4f}")
        logger.info(f"  Test Accuracy: {accuracy:.4f}")
        logger.info(f"  Test Precision: {precision:.4f}")
        logger.info(f"  Tuning time: {tuning_time:.2f} seconds")
        logger.info(f"  Best params: {best_params}")
        
        # Save individual model results to log
        with open('bacdive_tuning.log', 'a') as f:
            f.write(f"\n{model_name.upper()} Best Parameters:\n")
            for key, value in best_params.items():
                f.write(f"  {key}: {value}\n")
            f.write(f"CV score: {best_score:.4f}\n")
            f.write(f"Test F1: {f1:.4f}\n")
            f.write(f"Test Jaccard: {jaccard:.4f}\n")
            f.write(f"Test Accuracy: {accuracy:.4f}\n")
            f.write(f"Test Precision: {precision:.4f}\n")
            f.write("-"*40 + "\n")

    # Final summary
    logger.info("\n" + "="*80)
    logger.info("FINAL SUMMARY - ALL MODELS")
    logger.info("="*80)
    
    # Print results for all models
    for model_name in results['models'].keys():
        scores = results['best_scores'][model_name]
        params = results['best_params'][model_name]
        tuning_time = results['tuning_times'][model_name]
        
        logger.info(f"\n{model_name.upper()}:")
        logger.info(f"  Test F1 (macro): {scores['test_f1']:.4f}")
        logger.info(f"  Test Jaccard: {scores['test_jaccard']:.4f}")
        logger.info(f"  Test Accuracy: {scores['test_accuracy']:.4f}")
        logger.info(f"  Test Precision: {scores['test_precision']:.4f}")
        logger.info(f"  CV Jaccard: {scores['cv_score']:.4f}")
        logger.info(f"  Tuning time: {tuning_time:.2f}s")
        logger.info(f"  Best params: {params}")
    
    # Save results summary to file
    summary_path = os.path.join(save_dir, "tuning_summary.pkl")
    with open(summary_path, 'wb') as f:
        pickle.dump(results, f)
    logger.info(f"\nFull results summary saved to: {summary_path}")
    
    return results

# Usage: run on the feature-selected data
# results = build_and_tune_models(X_train_hybrid, y_train, X_test_hybrid, y_test, n_trials=30) 
results = build_and_tune_models(X_train_hybrid, y_train, X_test_hybrid, y_test, n_trials=120) # if only creating logistic regression models, can use a greater number of trials