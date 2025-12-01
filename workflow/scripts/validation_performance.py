import pandas as pd
import os
import sys
import joblib
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, accuracy_score
import re

# USAGE:
# python evaluate_models_validation.py <patho_dir> <benign_dir> <model_dir> <output_prefix>

pathodata = sys.argv[1]
benigndata = sys.argv[2]
model_dir = sys.argv[3]
output_prefix = sys.argv[4]

def normalize_chr_column(df):
    df = df.copy()
    df.iloc[:, 0] = df.iloc[:, 0].astype(str).str.replace("^chr", "", regex=True)
    df.iloc[:, 0] = df.iloc[:, 0].replace({"23": "X", "24": "Y"})
    return df

# Load data
pdel = pd.read_table(f"{pathodata}/{pathodata}bed_SBfinal100.bed", low_memory=False)
pneu = pd.read_table(f"{benigndata}/{benigndata}bed_SBfinal100.bed", low_memory=False)

pdel = normalize_chr_column(pdel)
pneu = normalize_chr_column(pneu)

# Label & features
pdel["copies"] = 1
pneu["copies"] = 1
pdel["label"] = 1
pneu["label"] = 0

# Match length
pdel["length"] = pdel.iloc[:, 2] - pdel.iloc[:, 1]
pneu["length"] = pneu.iloc[:, 2] - pneu.iloc[:, 1]

mean_len = pdel["length"].mean()
std_len = pdel["length"].std()
min_len, max_len = mean_len - std_len, mean_len + std_len

pdel_filtered = pdel[(pdel["length"] >= min_len) & (pdel["length"] <= max_len)]
pneu_filtered = pneu[(pneu["length"] >= min_len) & (pneu["length"] <= max_len)]

# Balance
if len(pdel_filtered) <= len(pneu_filtered):
    pdel_final = pdel_filtered
    pneu_final = pneu_filtered.sample(n=len(pdel_filtered), random_state=42)
else:
    pneu_final = pneu_filtered
    pdel_final = pdel_filtered.sample(n=len(pneu_filtered), random_state=42)

# Combined
all_sv = pd.concat([pdel_final, pneu_final]).reset_index(drop=True)
all_sv["key"] = all_sv.iloc[:, 0].astype(str) + "_" + all_sv.iloc[:, 1].astype(str) + "_" + all_sv.iloc[:, 2].astype(str)

X = all_sv.copy()
y = all_sv["label"].reset_index(drop=True)

# Pattern: only models with proper naming
model_pattern = re.compile(r"sunnytrainingGENjuly_(CB|XB)_\d+flank_cri.*\.joblib")

results = []

print("Evaluating models on validation set...")
print(os.listdir(model_dir))
for fname in os.listdir(model_dir):
    if not fname.endswith(".joblib") or not model_pattern.match(fname):
        continue

    model_path = os.path.join(model_dir, fname)
    print(model_path)
    try:
        model = joblib.load(model_path)
        features = model.feature_names_in_

        preds = model.predict_proba(X[features])[:, 1]
        auc = roc_auc_score(y, preds)
        acc = accuracy_score(y, preds > 0.5)

        results.append({
            "Model": fname,
            "AUC": auc,
            "Accuracy": acc
        })

        print(f"{fname} ? AUC: {auc:.4f} | Acc: {acc:.4f}")

    except Exception as e:
        print(f"Failed: {fname} - {e}")

# Save results
results_df = pd.DataFrame(results).sort_values(by="AUC", ascending=False)
csv_out = f"{output_prefix}_validation_results.csv"
results_df.to_csv(csv_out, index=False)
print(f"\nSaved to {csv_out}")

# Plot top 10
top = results_df.head(10)
plt.figure(figsize=(12, 6))
plt.barh(top["Model"], top["AUC"])
plt.xlabel("AUC")
plt.title("Top 10 Models on Validation Set (AUC)")
plt.gca().invert_yaxis()
plt.tight_layout()
plot_out = f"{output_prefix}_validation_plot.png"
plt.savefig(plot_out)
plt.show()
print(f"Plot saved to {plot_out}")
