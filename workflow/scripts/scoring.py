import os
import numpy as np
import pandas as pd
import joblib
import sys

from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    confusion_matrix,
    ConfusionMatrixDisplay,
    classification_report,
)

v2_SR_path = "workflow/models/cloudytrainingGENweighted_XB_100flank_crientropy_md30_mfsqrt_mln100_mid0.0_msl2_mss5_nes500.joblib"
v2_path = "workflow/models/sunnytrainingGENjuly_CB_100flank_crientropy_md30_mfsqrt_mln100_mid0.0_msl2_mss5_nes500.joblib"

rfSR = joblib.load(v2_SR_path)
rf = joblib.load(v2_path)


svdata = pd.read_table(sys.argv[2], low_memory=False)
copies_svdata = pd.read_table(sys.argv[1], header=None, low_memory=False)
svdata["copies"] = np.where(copies_svdata.iloc[:, 3].values == "DEL", -1, 1)


X_full = svdata[rfSR.feature_names_in_]

pred_v2 = rfSR.predict_proba(X_full)[:, 1]
pred_v2_lite = rf.predict_proba(X_full[rf.feature_names_in_])[:, 1]
copies_svdata.columns = ["chr", "start", "end", "type"]
copies_svdata["CADD-SV_score"] = pred_v2_lite
copies_svdata["CADD-SV-SR_score"] = pred_v2


output = pd.concat([copies_svdata, svdata.iloc[:,3:]], axis = 1)
output.to_csv(sys.argv[3], sep="\t", index=False)