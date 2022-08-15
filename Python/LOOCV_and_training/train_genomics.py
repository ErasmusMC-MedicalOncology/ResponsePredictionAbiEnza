import argparse
import pandas as pd
from sklearn.linear_model import LogisticRegression
import pickle

def init():
    parser = argparse.ArgumentParser(prog='train_genomics.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-input_file', help='input data matrix')
    parser.add_argument('-output_file', help='output')
    args = parser.parse_args()

    if not args.input_file or not args.output_file:
        parser.print_help()
        return 0
    return args

def run(args):
    train_cpct_df = pd.read_csv(args.input_file, index_col=[0], sep=",")
    train_cpct_df = train_cpct_df.drop(['Prior_AbiEnza', 'Prior_Chemo'], axis=1)

    # bad responder = 0, good responder = 1
    y_train_ = list(train_cpct_df['response'])
    y_train = [1 if "Good" in y else 0 for y in y_train_]

    train_cpct_df_input = train_cpct_df.drop('response', axis=1)

    X_train = train_cpct_df_input.values

    ### train logistic regression classifier
    logreg = LogisticRegression(solver='liblinear', random_state=1)
    logreg.fit(X_train, y_train)

    ### save model
    pickle.dump(logreg, open(args.output_file, 'wb'))

args = init()
if args is not 0:
    run(args)