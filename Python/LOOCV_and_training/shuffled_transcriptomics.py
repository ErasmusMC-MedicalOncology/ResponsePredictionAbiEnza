import sys
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneOut
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn import preprocessing


if __name__ == '__main__':
    input_folder = './raw/data/final_CPCT-02/'
    output_folder = './analysis/final_CPCT-02/classification_output/'
    cpct_df = pd.read_csv(str(input_folder + "training_set_transcriptomics_fastICA_40.csv"), index_col=[0], sep=",")

    y = list(cpct_df['response'])
    le = preprocessing.LabelEncoder()
    le.fit(y)
    y = le.transform(y)

    cpct_df_input = cpct_df.drop('response', axis=1)

    X = cpct_df_input.values

    ### load in permutated labels per round
    plabel_file = open(str(output_folder + "LOOCV_permuted_labels.csv"))
    permutated_lines = plabel_file.readlines()
    permutations = {}
    for px, pline in enumerate(permutated_lines):
        if px == 0:
            pheader = pline.strip("\n").split(",")
            test_samples = pheader[1:]
        else:
            pround = pline.strip("\n").split(",")
            permutations[pround[0]] = pround[1:]

    for px in permutations:
        p_list = permutations[px]
        shuffled_labels = [int(i) for i in p_list]

        y_perm = shuffled_labels
        y_perm = np.asarray(y_perm)
        y_perm = y_perm.astype(int)

        ### save selected features with test sample accuracy into a table
        output_file = open(
            str(output_folder + "LOOCV_linSVC_transcriptomics_40ICs_permutation_{).csv".format(
                px)), 'w')
        output_file.write("sampleId,true_label,prediction,class_0_probability,class_1_probability\n")

        ### initialize LOOCV
        loo = LeaveOneOut()
        n_splits = loo.get_n_splits(X)

        test_sample_index = 0
        test_samples = list(cpct_df_input.index)

        for train_index, test_index in loo.split(X):
            sample_id = test_samples[test_sample_index]
            test_sample_index += 1

            X_tv, X_test = X[train_index], X[test_index]
            y_tv, y_test = y_perm[train_index], y[test_index]

            ### train calibrated linSVC classifier with default settings
            l_svc = LinearSVC(C=1.0, penalty='l2', loss='squared_hinge', random_state=1, max_iter=10000)
            ccf = CalibratedClassifierCV(l_svc)

            ccf.fit(X_tv, y_tv)
            test_prob = ccf.predict_proba(X_test)
            pred_label = ccf.predict(X_test)
            test_p = test_prob.flatten()

            # predicted_label = 0 if class_0 > class_1 else 1
            output_file.write(
                "{},{},{},{},{}\n".format(sample_id, str(y_test[0]), str(pred_label[0]), str(test_p[0]),
                                          str(test_p[1])))