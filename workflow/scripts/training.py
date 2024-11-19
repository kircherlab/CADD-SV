import pandas as pd
import numpy as np
import os
import sys
from sklearn.ensemble import RandomForestClassifier   
import joblib
print("libraries loaded")
#NEED TO GENERATE COPIES AND DELETERIOUSNESS FOR EACH TS
input = pd.read_table(sys.argv[4], header=None, low_memory=False)
deleteriousness_feature = input[5]
copies_feature = [
    -1 if x == 'DEL' else 1 if x in ['DUP', 'INS'] else x 
    for x in input[3]
]
#print(deleteriousness_feature)


cb_data = pd.read_table(sys.argv[1], header=0, low_memory=False)
sb_data = pd.read_table(sys.argv[2], header=0, low_memory=False)
xb_data = pd.read_table(sys.argv[3], header=0, low_memory=False)

cb_data["copies"] = copies_feature
sb_data["copies"] = copies_feature
xb_data["copies"] = copies_feature

cb_data["deleteriousness"] = deleteriousness_feature
sb_data["deleteriousness"] = deleteriousness_feature
xb_data["deleteriousness"] = deleteriousness_feature
#LOAD 3 TS AND PUT THEM IN A LIST 
models = ["CB", "SB", "XB"]
TS_total = [cb_data, sb_data, xb_data]
n_estimators = [200]
min_samples_split = [5]
max_leaf_nodes = [100]
c = 0
for training in TS_total:
    print("starting " + models[c] + " training")
    y = training.pop("deleteriousness")
    if models[c] == "SB":
        X = training.iloc[:,4:]
    else:
        X = training.iloc[:,3:]
    X = pd.DataFrame(X, columns=X.columns, index=X.index)#fitted_scaler.transform(X), columns=X.columns, index=X.index)
    # Calculate the count of -1 and 1 in the 'copies' feature
    copy_1_count = (X['copies'] == 1).sum()
    copy_minus1_count = (X['copies'] == -1).sum()
    
    # Calculate the weight for copies == 1
    copy_1_weight = copy_minus1_count / copy_1_count
    
    # Assign sample weights based on the 'copies' feature
    sample_weights = np.where(X['copies'] == 1, copy_1_weight, 1)
    for estimator in n_estimators:
        for mss in min_samples_split:
            for mln in max_leaf_nodes:

                rf = RandomForestClassifier(n_estimators=estimator,
                        min_samples_split=mss, max_leaf_nodes=mln,
                        bootstrap=True, random_state=42)
                rf.fit(X, y)
                
                joblib.dump(rf, "./workflow/models/" + sys.argv[6] + models[c] + "_" + sys.argv[5] + "flank_" + ".joblib")
                print("./workflow/models/" + sys.argv[6] + models[c] + "_" + sys.argv[5] + "flank" + ".joblib")
                
                #joblib.dump(rf, "./workflow/models/" + sys.argv[6] + models[c] + "_" + sys.argv[5] + "flank_"
                #+ str(estimator) + "est_" + str(mss) + "_mss" +  str(mln)
                #+ "_mln" + models[c] + ".joblib")
                
    c += 1

