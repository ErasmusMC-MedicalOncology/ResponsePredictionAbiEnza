import sys
import os
import numpy as np
import argparse
import pandas as pd
import pickle

import scipy
from scipy.stats import ks_2samp

from precise.dimensionality_reduction import process_dim_reduction
from precise.principal_vectors import PVComputation
from precise.intermediate_factors import IntermediateFactors
from precise.pipeline_routine import ConsensusRepresentation

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import SparsePCA, PCA
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

def init():
    parser = argparse.ArgumentParser(prog='run_PRECISE_other.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-transcriptomics_validation_data_file', help='full file path to validation data input (see specified format in README.md)')
    parser.add_argument('-transcriptomics_training_data_file', help='full file path to "CPCT_cohort_transcriptomics_raw_data.csv"')
                          
    parser.add_argument('-output_folder', help='output folder name where prediction output and intermedier files are also saved')
    parser.add_argument('-prediction_output_file', help='output file name for prediction output')
    args = parser.parse_args()

    return args

def load_sparsePCA_components(training_transcriptomics_pc_model_file):
    training_transcriptomics_pc_model = load_model(training_transcriptomics_pc_model_file)

    t_loadings_ = training_transcriptomics_pc_model.components_
    t_sparsepc_loadings_ = scipy.linalg.orth(t_loadings_.transpose()).transpose()

    return t_sparsepc_loadings_

def run_sparsePCA(X_data):
    spca_model = SparsePCA(n_components=15, alpha=1)
    loadings_ = spca_model.fit(X_data).components_
    sparsepc_loadings_ = scipy.linalg.orth(loadings_.transpose()).transpose()

    return sparsepc_loadings_, spca_model

def standard_scale_data(data_to_scale):
    scaler = StandardScaler(with_mean=True, with_std=False)
    data_scaled = scaler.fit_transform(data_to_scale)
    scaled_df = pd.DataFrame(data_scaled, index=data_to_scale.index, columns=data_to_scale.columns)
    return scaled_df

def load_model(model_file_path):
    with open(model_file_path, 'rb') as file:
        model_ = pickle.load(file)
    return model_

def train_linSVC_classifier(X_training_transformed, y_training, output_folder):
    ### train calibrated linSVC classifier
    l_svc = LinearSVC(dual=False, random_state=1, max_iter=10000)
    ccf = CalibratedClassifierCV(l_svc)
    ccf.fit(X_training_transformed, y_training)

    # save model
    pickle.dump(ccf,
                open(os.path.join(output_folder, "trained_PRECISE_based_model.pkl"), 'wb'))

    return ccf

def run(args):
    training_set_with_response = pd.read_csv(args.transcriptomics_training_data_file, index_col=[0])
    validation_set_unsorted = pd.read_csv(args.transcriptomics_validation_data_file, index_col=[0])
    validation_samples = list(validation_set_unsorted.index.values)

    y_training = list(training_set_with_response['response'])
    training_set_unfiltered_unsorted = training_set_with_response.drop('response', axis=1)

    ### get gene list from validation file
    val_ENS_ids = set(list(validation_set_unsorted.columns.values))
    train_ENS_ids = set(list(training_set_unfiltered_unsorted.columns.values))
    common_ENS_ids = list(val_ENS_ids.intersection(train_ENS_ids))
    common_ENS_ids.sort()
    training_set = training_set_unfiltered_unsorted[common_ENS_ids]
    validation_set = validation_set_unsorted[common_ENS_ids]

    n_factors = 15
    n_pv = 15
    dim_reduction = 'sparsepca'
    dim_reduction_target = 'sparsepca'

    ### standard scale data
    X_training_data = training_set.values
    X_validation_data = validation_set.values

    ### consensus computation
    # use unscaled data as ConsensusRepresentation scales the data matrices
    consensus = ConsensusRepresentation(source_data=X_training_data,
                                        target_data=X_validation_data,
                                        n_factors=n_factors,
                                        n_pv=n_pv,
                                        dim_reduction='sparsepca',
                                        n_representations=15,
                                        use_data=False,
                                        mean_center=True,
                                        std_unit=False)

    consensus.fit(X_training_data)

    # save consensus model
    pickle.dump(consensus,
                open(os.path.join(args.output_folder, "consensus_model.pkl"), 'wb'))

    X_training_transformed = consensus.transform(X_training_data)
    X_validation_transformed = consensus.transform(X_validation_data)

    ### train classifier
    trained_model = train_linSVC_classifier(X_training_transformed, y_training, args.output_folder)

    ### predict validation cohort labels
    test_prob = trained_model.predict_proba(X_validation_transformed)
    pred_label = trained_model.predict(X_validation_transformed)
    test_p = test_prob

    # predicted_label = 0 if class_0 > class_1 else 1
    prediction_output_file = os.path.join(args.output_folder, args.prediction_output_file)
    prediction_output = open(prediction_output_file, 'w')
    prediction_output.write("sampleID,predicted_label,class_0_p,class_1_p\n")

    for vx, val_sample in enumerate(validation_samples):
        sample_test_p = test_p[vx]
        prediction_output.write("{},{},{},{}\n".format(val_sample, str(pred_label[vx]), str(sample_test_p[0]), str(sample_test_p[1])))

args = init()
if args is not 0:
    run(args)
