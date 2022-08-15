import argparse
import pandas as pd
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing

def init():
    parser = argparse.ArgumentParser(prog='LOOCV_only_chemo_arsi.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-input_file', help='input data matrix')
    parser.add_argument('-output_file', help='path to output with Leave-One-Out cross-validation results')
    args = parser.parse_args()

    if not args.input_file or not args.output_file:
        parser.print_help()
        return 0
    return args

def run(args):
    input_file = args.input_file
    output_file = args.output_file
    cpct_df = pd.read_csv(input_file, index_col=[0], sep=",")

    # bad responder = 0, good responder = 1
    y = list(cpct_df['response'])
    le = preprocessing.LabelEncoder()
    le.fit(y)
    y = le.transform(y)

    cpct_df_input = cpct_df.drop(['Genome.TMB','totalSV','SV.DUP','SV.DEL','response','Number_prior_treatment_lines'], axis=1)

    X = cpct_df_input.values

    ### save selected features with test sample accuracy into a table
    output_file = open(output_file, 'w')
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
        y_tv, y_test = y[train_index], y[test_index]

        ### train logistic regression classifier
        logreg = LogisticRegression(solver='liblinear', random_state=1)

        logreg.fit(X_tv, y_tv)
        test_prob = logreg.predict_proba(X_test)
        pred_label = logreg.predict(X_test)
        test_p = test_prob.flatten()

        # predicted_label = 0 if class_0 > class_1 else 1
        output_file.write(
            "{},{},{},{},{}\n".format(sample_id, str(y_test[0]), str(pred_label[0]), str(test_p[0]), str(test_p[1])))

args = init()
if args is not 0:
    run(args)