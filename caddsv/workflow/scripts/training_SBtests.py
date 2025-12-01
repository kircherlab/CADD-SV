import pandas as pd
import numpy as np
import os
import sys
import joblib
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import ParameterGrid

print("libraries loaded")

# Load input files
input = pd.read_table(sys.argv[4], header=None, low_memory=False)
deleteriousness_feature = input[5]
copies_feature = [-1 if x == 'DEL' else 1 if x in ['DUP', 'INS'] else x for x in input[3]]
chrom_feature = input[0].astype(str).str.replace("^chr", "", regex=True)

rsb_data = pd.read_table(sys.argv[1], header=0, low_memory=False)
sb_data = pd.read_table(sys.argv[2], header=0, low_memory=False)
db_data = pd.read_table(sys.argv[3], header=0, low_memory=False)
sbdb_data = pd.concat([sb_data, db_data.iloc[:,4:]], axis = 1)
#sbdbref_data = pd.concat([rsb_data, sbdb_data.iloc[:,3:]], axis = 1)

# Identify overlapping column names
overlap = rsb_data.columns.intersection(sbdb_data.columns)

# Create a copy of sbdb_data with renamed overlapping columns
sbdb_data_renamed = sbdb_data.rename(
    columns={col: f"{col}_alt" for col in overlap}
)

# Concatenate, taking only the desired columns from sbdb_data
sbdbref_data = pd.concat([rsb_data, sbdb_data_renamed.iloc[:, 4:]], axis=1)

# Add features
for data in [rsb_data, sb_data, db_data, sbdb_data, sbdbref_data]:
    data["copies"] = copies_feature
    data["deleteriousness"] = deleteriousness_feature
    data["chrom"] = chrom_feature

datasets = [rsb_data, sb_data, db_data, sbdb_data, sbdbref_data]
names = ["RSB", "SB", "DB", "SBDB", "RSBDB"]

# Hyperparameter grid (tuned for noisy data)
param_grid = {
    'n_estimators': [500],
    'min_samples_split': [5],
    'min_samples_leaf': [2],
    'max_leaf_nodes': [100],
    'max_depth': [30],
    'criterion': ['entropy'],
    'max_features': ['sqrt'],
    'min_impurity_decrease': [0.0]
}
results = []

os.makedirs("workflow/models", exist_ok=True)

for i, data in enumerate(datasets):
    filtered_data = data.copy()
    prefix = "GENweighted_"
    print(f"\nStarting {prefix}{names[i]} model training")

    train_data = filtered_data[filtered_data["chrom"] != "8"]
    holdout_data = filtered_data[filtered_data["chrom"] == "8"]

    y_train = train_data["deleteriousness"]
    y_holdout = holdout_data["deleteriousness"]

    if names[i] in ("RSB", "SB", "DB", "SBDB", "RSBDB"):
        X_train = train_data.iloc[:, 4:].copy()
        X_holdout = holdout_data.iloc[:, 4:].copy()
    else:
        X_train = train_data.iloc[:, 3:].copy()
        X_holdout = holdout_data.iloc[:, 3:].copy()

    for X in [X_train, X_holdout]:
        X.drop(columns=["deleteriousness", "chrom"], inplace=True, errors='ignore')
    copy_1_count = (X_train['copies'] == 1).sum()
    copy_minus1_count = (X_train['copies'] == -1).sum()
    
    # Calculate the weight for copies == 1
    copy_1_weight = copy_minus1_count / copy_1_count
    
    # Assign sample weights based on the 'copies' feature
    sample_weights = np.where(X_train['copies'] == 1, copy_1_weight, 1)
    for params in ParameterGrid(param_grid):
        rf = RandomForestClassifier(
            **params,
            bootstrap=True,
            random_state=42
        )
        rf.fit(X_train, y_train, sample_weight=sample_weights)

        if not X_holdout.empty:
            y_pred = rf.predict(X_holdout)
            acc = accuracy_score(y_holdout, y_pred)
        else:
            acc = None

        result_entry = {
            "Dataset": names[i],
            **params,
            "Holdout_Accuracy": acc
        }
        results.append(result_entry)

        # Custom abbreviation for each parameter
        abbreviations = {
            'n_estimators': 'nes',
            'min_samples_split': 'mss',
            'min_samples_leaf': 'msl',
            'max_leaf_nodes': 'mln',
            'max_depth': 'md',
            'criterion': 'cri',
            'max_features': 'mf',
            'min_impurity_decrease': 'mid'
        }

        param_str = "_".join([
            f"{abbreviations.get(k, k)}{v}" for k, v in params.items()
        ])

        model_filename = (
            f"{sys.argv[6]}{prefix}{names[i]}_{sys.argv[5]}flank_{param_str}.joblib"
        )
        model_path = os.path.join("workflow/models", model_filename)
        joblib.dump(rf, model_path)
        print(f"Saved model: {model_filename} | Holdout Acc: {acc:.4f}" if acc is not None else "No holdout data")

# Convert results to DataFrame
results_df = pd.DataFrame(results)
csv_path = "workflow/models/holdout_results.csv"
results_df.to_csv(csv_path, index=False)
print(f"\nAll results saved to {csv_path}")

# Plot top 10 models by holdout accuracy
top_df = results_df.dropna(subset=["Holdout_Accuracy"]).sort_values(by="Holdout_Accuracy", ascending=False).head(10)

plt.figure(figsize=(12, 6))
plt.barh(
    [f"{row.Dataset}_{i+1}" for i, row in top_df.iterrows()],
    top_df["Holdout_Accuracy"]
)
plt.xlabel("Holdout Accuracy")
plt.title("Top 10 Models by Holdout Accuracy (chr8)")
plt.gca().invert_yaxis()
plt.tight_layout()
plot_path = "workflow/models/top_models_accuracy_plot.png"
plt.savefig(plot_path)
plt.show()
print(f"Accuracy plot saved to {plot_path}")
