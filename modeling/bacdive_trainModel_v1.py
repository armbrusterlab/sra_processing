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
terms_keep = sorted([k for k in terms.keys() if terms[k] >= minSize])

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
# Initialize TF-IDF Vectorizer
vectorizer = TfidfVectorizer(stop_words='english', max_df=0.7) # first-pass "feature selection"; remove words that appear in >70% of the records since they're not helpful
# I could use min_df but opted not to, as many groups are very small and this might remove legitimate information

# Transform the text data to feature vectors
# first, convert the ### delimiter in the isolation source strings into spaces so that the words on either side of the delimiter will be processed properly
# also remove the <I> and </I> for italicization
iso_source = [
    row["Isolation source"].replace("###", " ").replace("<I>", "").replace("</I>", "")
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
selected_indices_hybrid = hybrid_feature_selection(
    X_train, y_train, chi2_k=4000, rf_k=2000
)
X_train_hybrid = X_train[:, selected_indices_hybrid]
X_test_hybrid = X_test[:, selected_indices_hybrid]
# These dataframes contain the top 2000 features.
# >>> X_train_hybrid.shape  
# (34626, 2000)
# >>> X_test_hybrid.shape 
# (8657, 2000)

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
from sklearn.metrics import f1_score, jaccard_score, accuracy_score
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
    
    # Dictionary to store the best model overall
    best_overall = {
        'model_name': None, 
        'jaccard': -1, 
        'model': None,
        'f1': None,
        'accuracy': None,
        'params': None
    }
    
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
        }#,
        # 'stacking': {
        #     'final_estimator__C': optuna.distributions.FloatDistribution(0.1, 10.0),
        #     'final_estimator__penalty': optuna.distributions.CategoricalDistribution(['l1', 'l2']),
        #     'final_estimator__solver': optuna.distributions.CategoricalDistribution(['liblinear']),
        #     'final_estimator__max_iter': optuna.distributions.CategoricalDistribution([1000])
        # }
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
        lr = LogisticRegression(
            penalty=params['estimator__penalty'],
            C=params['estimator__C'],
            solver=params['estimator__solver'],
            max_iter=params['estimator__max_iter'],
            random_state=42
        )
        return MultiOutputClassifier(lr, n_jobs=-1)
    
    # def create_stacking(params, base_models):
    #     """Create stacking classifier with pre-trained base models."""
    #     meta_lr = LogisticRegression(
    #         C=params['final_estimator__C'],
    #         penalty=params['final_estimator__penalty'],
    #         solver=params['final_estimator__solver'],
    #         max_iter=params['final_estimator__max_iter'],
    #         random_state=42
    #     )
        
    #     stack = StackingClassifier(
    #         estimators=base_models,
    #         final_estimator=meta_lr,
    #         cv=5,
    #         stack_method='predict_proba',
    #         n_jobs=-1
    #     )
    #     return stack
    
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
            model = create_model_func(params)
            
            # Cross-validation with Jaccard score
            cv = KFold(n_splits=5, shuffle=True, random_state=42)
            
            # Compute cross-validation scores
            scores = cross_val_score(
                model, X_train, y_train,
                cv=cv,
                scoring=jaccard_scorer,
                n_jobs=-1,
                error_score='raise'
            )
            
            return scores.mean()
            
        except Exception as e:
            logger.warning(f"Trial failed: {str(e)}")
            return -1.0
    
    def save_best_model(model, model_name, params, scores, save_dir):
        """Save the best performing model to disk."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save the best model
        model_path = os.path.join(save_dir, "best_model.joblib")
        joblib.dump(model, model_path, compress=3)
        
        # Save metadata
        metadata = {
            'model_name': model_name,
            'timestamp': timestamp,
            'parameters': params,
            'scores': scores,
            'feature_count': X_train.shape[1],
            'label_count': y_train.shape[1],
            'training_date': datetime.now().isoformat()
        }
        
        metadata_path = os.path.join(save_dir, "best_model_metadata.pkl")
        with open(metadata_path, 'wb') as f:
            pickle.dump(metadata, f)
        
        logger.info(f"\n? Best model saved to: {model_path}")
        logger.info(f"? Metadata saved to: {metadata_path}")
        
        return model_path, metadata_path
    
    # Optimize each model type
    model_configs = [
        ('random_forest', create_random_forest, param_distributions['random_forest']),
        ('gradient_boosting', create_gradient_boosting, param_distributions['gradient_boosting']),
        ('lasso', create_lasso, param_distributions['lasso']),
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
        
        scores = {
            'cv_score': best_score,
            'test_f1': f1,
            'test_jaccard': jaccard,
            'test_accuracy': accuracy
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
        
        # Update best overall if this model is better
        if jaccard > best_overall['jaccard']:
            best_overall['jaccard'] = jaccard
            best_overall['model_name'] = model_name
            best_overall['model'] = best_model
            best_overall['f1'] = f1
            best_overall['accuracy'] = accuracy
            best_overall['params'] = best_params
            best_overall['cv_score'] = best_score
        
        # Print results immediately
        logger.info(f"\n{model_name.upper()} Results:")
        logger.info(f"  Best CV Jaccard: {best_score:.4f}")
        logger.info(f"  Test F1 (macro): {f1:.4f}")
        logger.info(f"  Test Jaccard: {jaccard:.4f}")
        logger.info(f"  Test Accuracy: {accuracy:.4f}")
        logger.info(f"  Tuning time: {tuning_time:.2f} seconds")
        logger.info(f"  Best params: {best_params}")
        
        # Save individual model results to log
        with open('bacdive_tuning.log', 'a') as f:
            f.write(f"\n{model_name.upper()} Best Parameters:\n")
            for key, value in best_params.items():
                f.write(f"  {key}: {value}\n")
            f.write(f"CV Jaccard: {best_score:.4f}\n")
            f.write(f"Test F1: {f1:.4f}\n")
            f.write(f"Test Jaccard: {jaccard:.4f}\n")
            f.write(f"Test Accuracy: {accuracy:.4f}\n")
            f.write("-"*40 + "\n")
    
    # # Now optimize Stacking classifier
    # logger.info(f"\n{'='*60}")
    # logger.info("Optimizing STACKING CLASSIFIER")
    # logger.info(f"{'='*60}")
    
    # # Retrain best base models if they are MultiOutputClassifier
    # stacking_base_models = []
    # for name, model in results['models'].items():
    #     if name != 'lasso' and model is not None:
    #         if isinstance(model, MultiOutputClassifier):
    #             stacking_base_models.append((name, clone(model.estimator)))
    #         else:
    #             stacking_base_models.append((name, clone(model)))
    
    # # Optimize stacking
    # start_time = time.time()
    
    # def stacking_objective(trial):
    #     """Objective function for stacking optimization."""
    #     params = {}
    #     for param_name, dist in param_distributions['stacking'].items():
    #         if isinstance(dist, optuna.distributions.CategoricalDistribution):
    #             params[param_name] = trial.suggest_categorical(param_name, dist.choices)
    #         elif isinstance(dist, optuna.distributions.IntDistribution):
    #             params[param_name] = trial.suggest_int(param_name, dist.low, dist.high)
    #         elif isinstance(dist, optuna.distributions.FloatDistribution):
    #             params[param_name] = trial.suggest_float(param_name, dist.low, dist.high, log=dist.log)
        
    #     try:
    #         # Create stacking model
    #         meta_lr = LogisticRegression(
    #             C=params['final_estimator__C'],
    #             penalty=params['final_estimator__penalty'],
    #             solver=params['final_estimator__solver'],
    #             max_iter=params['final_estimator__max_iter'],
    #             random_state=42
    #         )
            
    #         # Ensure base models are cloned
    #         base_clones = []
    #         for name, estimator in stacking_base_models:
    #             base_clones.append((name, clone(estimator)))
            
    #         stack = StackingClassifier(
    #             estimators=base_clones,
    #             final_estimator=meta_lr,
    #             cv=5,
    #             stack_method='predict_proba',
    #             n_jobs=-1
    #         )
            
    #         # Cross-validation for stacking
    #         cv = KFold(n_splits=5, shuffle=True, random_state=42)
    #         scores = []
    #         for train_idx, val_idx in cv.split(X_train):
    #             X_train_fold, X_val_fold = X_train[train_idx], X_train[val_idx]
    #             y_train_fold, y_val_fold = y_train[train_idx], y_train[val_idx]
                
    #             # Fit stacking on fold
    #             stack_clone = clone(stack)
    #             stack_clone.fit(X_train_fold, y_train_fold)
    #             y_pred_fold = stack_clone.predict(X_val_fold)
    #             score = jaccard_score(y_val_fold, y_pred_fold, average='samples', zero_division=0)
    #             scores.append(score)
            
    #         return np.mean(scores)
            
    #     except Exception as e:
    #         logger.warning(f"Stacking trial failed: {str(e)}")
    #         return -1.0
    
    # stacking_study = optuna.create_study(
    #     direction='maximize',
    #     sampler=optuna.samplers.RandomSampler(seed=42)
    # )
    
    # stacking_study.optimize(
    #     stacking_objective,
    #     n_trials=n_trials,
    #     show_progress_bar=True
    # )
    
    # tuning_time = time.time() - start_time
    
    # # Train best stacking model
    # best_stack_params = stacking_study.best_params
    # best_stack_score = stacking_study.best_value
    
    # logger.info(f"\nTraining best stacking model on full training data...")
    
    # # Create final stacking model with best params
    # meta_lr = LogisticRegression(
    #     C=best_stack_params['final_estimator__C'],
    #     penalty=best_stack_params['final_estimator__penalty'],
    #     solver=best_stack_params['final_estimator__solver'],
    #     max_iter=best_stack_params['final_estimator__max_iter'],
    #     random_state=42
    # )
    
    # # Retrain base estimators on full data
    # base_estimators_full = []
    # for name, estimator in stacking_base_models:
    #     base_estimators_full.append((name, clone(estimator)))
    
    # best_stacking = MultiOutputClassifier(
    #     StackingClassifier(
    #         estimators=base_estimators_full,
    #         final_estimator=meta_lr,
    #         cv=5,
    #         stack_method='predict_proba',
    #         n_jobs=-1
    #     ),
    #     n_jobs=-1
    # )
    # best_stacking.fit(X_train, y_train)
    
    # # Evaluate stacking
    # y_pred_stack = best_stacking.predict(X_test)
    # f1_stack = f1_score(y_test, y_pred_stack, average='macro', zero_division=0)
    # jaccard_stack = jaccard_score(y_test, y_pred_stack, average='samples', zero_division=0)
    # accuracy_stack = accuracy_score(y_test, y_pred_stack)
    
    # scores_stack = {
    #     'cv_score': best_stack_score,
    #     'test_f1': f1_stack,
    #     'test_jaccard': jaccard_stack,
    #     'test_accuracy': accuracy_stack
    # }
    
    # # Store stacking results
    # results['models']['stacking'] = best_stacking
    # results['best_params']['stacking'] = best_stack_params
    # results['best_scores']['stacking'] = scores_stack
    # results['tuning_times']['stacking'] = tuning_time
    
    # # Update best overall if stacking is better
    # if jaccard_stack > best_overall['jaccard']:
    #     best_overall['jaccard'] = jaccard_stack
    #     best_overall['model_name'] = 'stacking'
    #     best_overall['model'] = best_stacking
    #     best_overall['f1'] = f1_stack
    #     best_overall['accuracy'] = accuracy_stack
    #     best_overall['params'] = best_stack_params
    #     best_overall['cv_score'] = best_stack_score
    
    # logger.info(f"\nSTACKING Results:")
    # logger.info(f"  Best CV Jaccard: {best_stack_score:.4f}")
    # logger.info(f"  Test F1 (macro): {f1_stack:.4f}")
    # logger.info(f"  Test Jaccard: {jaccard_stack:.4f}")
    # logger.info(f"  Test Accuracy: {accuracy_stack:.4f}")
    # logger.info(f"  Tuning time: {tuning_time:.2f} seconds")
    # logger.info(f"  Best params: {best_stack_params}")
    
    # # Save stacking results to log
    # with open('bacdive_tuning.log', 'a') as f:
    #     f.write(f"\nSTACKING Best Parameters:\n")
    #     for key, value in best_stack_params.items():
    #         f.write(f"  {key}: {value}\n")
    #     f.write(f"CV Jaccard: {best_stack_score:.4f}\n")
    #     f.write(f"Test F1: {f1_stack:.4f}\n")
    #     f.write(f"Test Jaccard: {jaccard_stack:.4f}\n")
    #     f.write(f"Test Accuracy: {accuracy_stack:.4f}\n")
    #     f.write("="*40 + "\n")
    
    # Save the best performing model
    logger.info(f"\n{'='*60}")
    logger.info("Saving Best Performing Model")
    logger.info(f"{'='*60}")
    
    model_path, metadata_path = save_best_model(
        best_overall['model'],
        best_overall['model_name'],
        best_overall['params'],
        {
            'cv_score': best_overall.get('cv_score', 0),
            'test_f1': best_overall['f1'],
            'test_jaccard': best_overall['jaccard'],
            'test_accuracy': best_overall['accuracy']
        },
        save_dir
    )
    
    # Update results with best overall
    results['best_overall'] = best_overall
    results['best_overall']['saved_path'] = model_path
    results['best_overall']['saved_metadata_path'] = metadata_path
    
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
        logger.info(f"  CV Jaccard: {scores['cv_score']:.4f}")
        logger.info(f"  Tuning time: {tuning_time:.2f}s")
        logger.info(f"  Best params: {params}")
    
    # Best overall model
    logger.info("\n" + "="*80)
    logger.info(f"?? BEST OVERALL MODEL: {best_overall['model_name'].upper()}")
    logger.info(f"  Test F1 (macro): {best_overall['f1']:.4f}")
    logger.info(f"  Test Jaccard: {best_overall['jaccard']:.4f}")
    logger.info(f"  Test Accuracy: {best_overall['accuracy']:.4f}")
    logger.info(f"  CV Jaccard: {best_overall.get('cv_score', 0):.4f}")
    logger.info(f"  Saved to: {model_path}")
    logger.info("="*80)
    
    # Save results summary to file
    summary_path = os.path.join(save_dir, "tuning_summary.pkl")
    with open(summary_path, 'wb') as f:
        pickle.dump(results, f)
    logger.info(f"\nFull results summary saved to: {summary_path}")
    
    return results

# Helper function to load the saved best model
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

# # Usage example:
# if __name__ == "__main__":
#     # Assuming X_train, X_test, y_train, y_test are already defined
    
#     # Run tuning and save only the best model
#     results = build_and_tune_models(
#         X_train, y_train, X_test, y_test, 
#         n_trials=10, 
#         save_dir='models'
#     )
    
#     # Later, load the best model:
#     model, metadata = load_best_model('models')
#     print(f"Loaded model: {metadata['model_name']}")
#     print(f"Jaccard score: {metadata['scores']['test_jaccard']:.4f}")
#     print(f"F1 score: {metadata['scores']['test_f1']:.4f}")
#     print(f"Parameters: {metadata['parameters']}")
    
#     # Make predictions with loaded model
#     y_pred = model.predict(X_test)

# Usage: run on the feature-selected data
results = build_and_tune_models(X_train_hybrid, y_train, X_test_hybrid, y_test, n_trials=30) # for the real thing: use the feature-selected X

# # TODO start: to test out the following code I will subsample rows from the training dataset, but for the real thing I will need to comment out the following code.
# # also arbitrarily cut down on the number of columns so as to expedite the testing
# train_idx_subsample = train_idx[:1000]
# X_train = X[train_idx_subsample, :100]
# y_train = y.iloc[train_idx_subsample]
# X_test = X_test[:, :100]
# # Also note that this was based on the original X, but for the real thing I'll use the feature-selected X_train_hybrid and X_test_hybrid.
# results = build_and_tune_models(X_train, y_train, X_test, y_test, n_trials=10)
# # turns out that by subsampling the data, I ended up with some classes that had <2 positive samples, i.e. some of the 65 category-subcategory pairs didn't appear in the subsample.
# # This caused errors with gradient boosting, but it shouldn't occur as long as I don't subsample the data.
# # TODO end

