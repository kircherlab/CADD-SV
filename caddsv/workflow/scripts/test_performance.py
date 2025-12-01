import pandas as pd
import os
import joblib
import sys
import re
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score

# USAGE:
# python evaluate_models_chr8.py cb_data.txt xb_data.txt input.txt output_prefix

cb_path = sys.argv[1]
xb_path = sys.argv[2]
input_path = sys.argv[3]
output_prefix = sys.argv[4]

# Load input and holdout data
input_df = pd.read_table(input_path, header=None, low_memory=False)
deleteriousness_feature = input_df[5]
copies_feature = [-1 if x == 'DEL' else 1 if x in ['DUP', 'INS'] else x for x in input_df[3]]
chrom_feature = input_df[0].astype(str).str.replace("^chr", "", regex=True)

cb_data = pd.read_table(cb_path, header=0, low_memory=False)
xb_data = pd.read_table(xb_path, header=0, low_memory=False)

for df in [cb_data, xb_data]:
    df["copies"] = copies_feature
    df["deleteriousness"] = deleteriousness_feature
    df["chrom"] = chrom_feature

datasets = {
    "CB": cb_data,
    "XB": xb_data
}

results = []
model_dir = "workflow/models"

# Regex pattern for models matching the format
pattern = re.compile(r"GENjuly_(CB|XB)_\d+flank_cri.*\.joblib")

print("Evaluating models on chr8 holdout...")

for filename in os.listdir(model_dir):
    if filename.endswith(".joblib") and pattern.search(filename):
        dataset_key = "CB" if "_CB_" in filename else "XB" if "_XB_" in filename else None
        if dataset_key is None:
            continue

        df = datasets[dataset_key]
        holdout = df[df["chrom"] == "8"].copy()
        y_holdout = holdout["deleteriousness"]

        if dataset_key == "CB":
            X_holdout = holdout.iloc[:, 3:].copy()
        else:  # XB
            X_holdout = holdout.iloc[:, 3:].copy()

        X_holdout.drop(columns=["deleteriousness", "copies", "chrom"], inplace=True, errors='ignore')

        model_path = os.path.join(model_dir, filename)
        try:
            model = joblib.load(model_path)
            y_pred = model.predict(X_holdout)
            acc = accuracy_score(y_holdout, y_pred)

            results.append({
                "Dataset": dataset_key,
                "Model": filename,
                "Accuracy": acc
            })

            print(f"{filename}: {acc:.4f}")

        except Exception as e:
            print(f"Error loading or evaluating {filename}: {e}")

# Save results
results_df = pd.DataFrame(results)
csv_out = f"{output_prefix}_chr8_holdout_accuracy.csv"
results_df.to_csv(csv_out, index=False)
print(f"\nSaved evaluation results to {csv_out}")

# Plot top 10 models
top_df = results_df.sort_values(by="Accuracy", ascending=False).head(10)

plt.figure(figsize=(12, 6))
plt.barh(top_df["Model"], top_df["Accuracy"])
plt.xlabel("Accuracy on chr8")
plt.title("Top 10 Models on chr8 Holdout")
plt.gca().invert_yaxis()
plt.tight_layout()
plot_out = f"{output_prefix}_chr8_top10_plot.png"
plt.savefig(plot_out)
plt.show()
print(f"Saved plot to {plot_out}")
