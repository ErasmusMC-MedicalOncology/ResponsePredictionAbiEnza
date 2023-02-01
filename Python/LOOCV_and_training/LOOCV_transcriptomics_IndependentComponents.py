import argparse
import pandas as pd
from sklearn.model_selection import LeaveOneOut
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn import preprocessing

def init():
    parser = argparse.ArgumentParser(prog='LOOCV_transcriptomics.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-transcriptomics_input_folder', help='folder with transcriptomics matrix')
    parser.add_argument('-output_file', help='path to output with Leave-One-Out cross-validation results')
    parser.add_argument('-component', help='number of independent components (ICs) to use')
    args = parser.parse_args()

    if not args.genomics_input_file or not args.transcriptomics_input_file or not args.component or not args.output_file:
        parser.print_help()
        return 0
    return args

def run(args):
    component = args.component
    transcriptomics_input_folder = args.transcriptomics_input_folder
    output_file = args.output_file
    cpct_df = pd.read_csv(str(transcriptomics_input_folder + "training_set_transcriptomics_ICs_{}.csv".format(component)), index_col=[0], sep=",")

    y = list(cpct_df['response'])
    le = preprocessing.LabelEncoder()
    le.fit(y)
    y = le.transform(y)

    cpct_df_input = cpct_df.drop('response', axis=1)

    X = cpct_df_input.values

    output_file = open(output_file, 'w')
    output_file.write("sampleId,true_label,prediction,class_0_probability,class_1_probability\n")

    ### initialize LOOCV
    loo = LeaveOneOut()

    test_sample_index = 0
    test_samples = list(cpct_df_input.index)

    for train_index, test_index in loo.split(X):
        sample_id = test_samples[test_sample_index]
        test_sample_index += 1

        X_tv, X_test = X[train_index], X[test_index]
        y_tv, y_test = y[train_index], y[test_index]

        l_svc = LinearSVC(dual=False, random_state=1, max_iter=10000)

        ccf = CalibratedClassifierCV(l_svc)

        ccf.fit(X_tv, y_tv)
        test_prob = ccf.predict_proba(X_test)
        pred_label = ccf.predict(X_test)
        test_p = test_prob.flatten()

        # predicted_label = 0 if class_0 > class_1 else 1
        output_file.write(
            "{},{},{},{},{}\n".format(sample_id, str(y_test[0]), str(pred_label[0]), str(test_p[0]), str(test_p[1])))

args = init()
if args is not 0:
    run(args)
