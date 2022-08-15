import sys
import argparse
import pandas as pd
from sklearn.model_selection import LeaveOneOut

from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

from sklearn.preprocessing import StandardScaler
from sklearn import preprocessing

from sklearn.ensemble import StackingClassifier

def init():
    parser = argparse.ArgumentParser(prog='LOOCV_stacking.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-genomics_input_file', help='input genomics matrix')
    parser.add_argument('-transcriptomics_input_folder', help='folder with transcriptomics matrix')
    parser.add_argument('-output_file', help='path to output with Leave-One-Out cross-validation results')
    parser.add_argument('-component', help='number of sPCA components to use')
    args = parser.parse_args()

    if not args.genomics_input_file or not args.transcriptomics_input_file or not args.component or not args.output_file:
        parser.print_help()
        return 0
    return args

def run(args):
    component = args.component
    genomics_input_file = args.genomics_input_file
    transcriptomics_input_folder = args.transcriptomics_input_folder
    output_file = args.output_file
    transcriptomics_df = pd.read_csv(str(transcriptomics_input_folder + "training_set_transcriptomics_sparsePCA_{}.csv".format(component)), index_col=[0], sep=",")
    genomics_df = pd.read_csv(genomics_input_file, index_col=[0], sep=",")

    # bad responder = 0, good responder = 1
    y = list(genomics_df['response'])
    le = preprocessing.LabelEncoder()
    le.fit(y)
    y = le.transform(y)

    genomics_df_ = genomics_df.drop(['response', 'Prior_Chemo', 'Prior_AbiEnza', 'Number_prior_treatment_lines'], axis=1)
    transcriptomics_df_ = transcriptomics_df.drop('response', axis=1)

    merged_input = pd.merge(genomics_df_, transcriptomics_df_, left_index=True, right_index=True)

    X = merged_input.values

    #### Genomics classifier first
    output_file = open(output_file, 'w')
    output_file.write("sampleId,true_label,prediction,class_0_p,class_1_p\n")

    ### initialize LOOCV
    loo = LeaveOneOut()

    test_sample_index = 0
    test_samples = list(merged_input.index)

    # center data
    scaler = StandardScaler(with_mean=True, with_std=False)
    X_scaled = scaler.fit_transform(X)

    for train_index, test_index in loo.split(X_scaled):
        sample_id = test_samples[test_sample_index]
        test_sample_index += 1

        X_tv, X_test = X_scaled[train_index], X_scaled[test_index]
        y_tv, y_test = y[train_index], y[test_index]

        estimators = [('svr', CalibratedClassifierCV(LinearSVC(random_state=1, max_iter=10000))),
                      ('logreg', LogisticRegression(solver='liblinear', random_state=1))]

        clf = StackingClassifier(estimators = estimators, final_estimator = LogisticRegression())

        clf.fit(X_tv, y_tv)
        test_prob = clf.predict_proba(X_test)
        pred_label = clf.predict(X_test)
        test_p = test_prob.flatten()

        output_file.write(
            "{},{},{},{},{}\n".format(sample_id, str(y_test[0]), str(pred_label[0]), str(test_p[0]), str(test_p[1])))

args = init()
if args is not 0:
    run(args)