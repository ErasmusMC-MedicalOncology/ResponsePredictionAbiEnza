import sys
import argparse
import pandas as pd
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer

def init():
    parser = argparse.ArgumentParser(prog='genomics_w_clinvar_prediction.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-genomics_clinvar_validation_data_file', help='full file path to validation data input (see specified format in README.md)')
    parser.add_argument('-genomics_model', help='full path of "genomics.pkl"')
    parser.add_argument('-genomics_clinvar_model', help='full path of "genomics_w_clinvar.pkl"')
    parser.add_argument('-output', help='full path of output file where raw predictions are saved')
    args = parser.parse_args()

    if not args.genomics_clinvar_validation_data_file or not args.genomics_model \
            or not args.genomics_clinvar_model or not args.output:
        parser.print_help()
        return 0
    return args

def load_model(model_file_path):
    with open(model_file_path, 'rb') as file:
        model_ = pickle.load(file)
    return model_

def run(args):
    genomics_clinvar_validation_data_file = args.genomics_clinvar_validation_data_file
    genomics_model_file = args.genomics_model
    genomics_clinvar_model_file = args.genomics_clinvar_model
    prediction_outfile = args.output

    ### features
    genomic_features = ['Genome.TMB', 'totalSV', 'SV.DUP', 'SV.DEL']
    clinvar_features = ['Prior_AbiEnza']

    df = pd.read_csv(genomics_clinvar_validation_data_file, index_col=[0])
    #dm = df.T

    ### filter input data
    X_gen = df[genomic_features]
    X_clinvar = df[clinvar_features]

    ### center data
    #scaler_gen = StandardScaler(with_mean=True, with_std=True)
    #X_gen_scaled = pd.DataFrame(scaler_gen.fit_transform(X_gen), columns=X_gen.columns, index=X_gen.index)
    X_gen_clinvar_scaled = X_gen.join(X_clinvar)

    ### Load models
    genomics_model = load_model(genomics_model_file)
    genomics_clinvar_model = load_model(genomics_clinvar_model_file)

    ### Make predictions
    gen_prob = genomics_model.predict_proba(X_gen)
    gen_labels = genomics_model.predict(X_gen)

    gen_clinvar_prob = genomics_clinvar_model.predict_proba(X_gen_clinvar_scaled)
    gen_clinvar_labels = genomics_clinvar_model.predict(X_gen_clinvar_scaled)

    ### save predictions
    prediction_output = open(prediction_outfile, 'w')
    prediction_output.write(
        "sampleID,true_label,pred_label_genomics,p_genomics_class_0,p_genomics_class_1,pred_label_genomics_w_clinvar,"
        "p_genomics_w_clinvar_class_0,p_genomics_w_clinvar_class_1\n")
    sample_ids = [s for s in list(df.index)]
    labels = list(df['response'])
    for sx, sample_ in enumerate(sample_ids):
        label = labels[sx]
        p_gen = gen_prob[sx]
        label_gen = gen_labels[sx]

        p_gen_clinvar = gen_clinvar_prob[sx]
        label_gen_clinvar = gen_clinvar_labels[sx]
        prediction_output.write(
            "{},{},{},{},{},{},{},{}\n".format(sample_, label, label_gen, p_gen[0], p_gen[1], label_gen_clinvar, p_gen_clinvar[0], p_gen_clinvar[1]))

args = init()
if args is not 0:
    run(args)
