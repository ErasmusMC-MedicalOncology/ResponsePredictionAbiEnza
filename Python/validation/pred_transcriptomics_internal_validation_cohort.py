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
    input_folder = '/hpc/compgen/projects/abi-enza_response_prediction/response_group_classification/raw/data/final_CPCT-02/'
    output_folder = '/hpc/compgen/projects/abi-enza_response_prediction/response_group_classification/analysis/adanyi/final_CPCT-02/classification_output/'
    cpct_df = pd.read_csv(str(input_folder + "internal_cohort_transcriptomics_sparsePCA_15.csv"), index_col=[0], sep=",")
    training_cpct_df = pd.read_csv(str(input_folder + "training_set_transcriptomics_sparsePCA_15.csv"), index_col=[0], sep=",")

    training_cpct_df = training_cpct_df[training_cpct_df['response'] != "Ambiguous Responder (101-179 days)"]

    y = list(training_cpct_df['response'])
    le = preprocessing.LabelEncoder()
    le.fit(y)
    y_trans = le.transform(y)

    cpct_df_input = cpct_df.drop('response', axis=1)
    training_cpct_df_input = training_cpct_df.drop('response', axis=1)

    y_test = list(cpct_df['response'])

    X_train = training_cpct_df_input.values
    X_test = cpct_df_input.values

    ### save selected features with test sample accuracy into a table
    output_file = open(str(output_folder + "pred_internal_transcriptomics_sparsePCA_15.csv"), 'w')
    output_file.write("sampleId,true_label,prediction,class_0_probability,class_1_probability\n")

    l_svc = LinearSVC(dual=False, random_state=1, max_iter=10000)
    ccf = CalibratedClassifierCV(l_svc)

    ccf.fit(X_train, y_trans)
    test_prob = ccf.predict_proba(X_test)
    pred_label = ccf.predict(X_test)

    sample_ids = [s for s in list(cpct_df_input.index)]
    for sx, sample_ in enumerate(sample_ids):
        label = y_test[sx]
        p_gen = test_prob[sx]
        label_gen = pred_label[sx]
        output_file.write(
            "{},{},{},{},{}\n".format(sample_, label, label_gen, p_gen[0], p_gen[1]))
