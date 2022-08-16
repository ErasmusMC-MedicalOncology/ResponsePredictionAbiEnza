import os
import pandas as pd
import argparse
import pickle

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import SparsePCA, PCA

import sys

def init():
    parser = argparse.ArgumentParser(prog='train_genomics_arsi.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-metadata_path', help='metadata file path')
    parser.add_argument('-raw_data', help='raw data path')
    parser.add_argument('-output_dir', help='output directory')
    parser.add_argument('-component', help='number of components')
    args = parser.parse_args()

    if not args.metadata_path or not args.raw_data or not args.output_dir or not args.component:
        parser.print_help()
        return 0
    return args

def transform_with_fitted_sparsePCA(model_path, data):
    base_path = model_path.split(".")[0].split("_")
    component_number = int(base_path[len(base_path)-1])
    model = pickle.load(open(model_path, "rb"))
    data_without_response = data.drop("response", axis=1)
    t_data = model.transform(data_without_response.values)
    transformed_df = pd.DataFrame(t_data, index=data.index,
                                  columns=["C{}".format(comp_) for comp_ in range(1, component_number+1)])
    transformed_df['response'] = list(data['response'])

    return transformed_df, component_number


def run_sparsePCA(data_type, data_set, components, output_dir):
    n_factors = components

    data_set_unlabelled = data_set.drop('response', axis=1)
    X_data = data_set_unlabelled.values

    # Compute source PCs
    pc_model = SparsePCA(n_factors, alpha=1)
    pc_model.fit(X_data)
    X_transformed = pc_model.transform(X_data)

    # save model
    pickle.dump(pc_model,
                open(os.path.join(output_dir, '{}_sparsePCA_{}.pkl'.format(data_type, components), 'wb')))

    # save transformed data matrix
    transformed_df = pd.DataFrame(X_transformed, index=data_set_unlabelled.index,
                                  columns=["C{}".format(comp_) for comp_ in range(1, components+1)])
    transformed_df.to_csv(
        os.path.join(output_dir, "sparsePCA_components_{}_{}.csv".format(data_type, components, header=True, index=True)))

    return transformed_df


def standard_scale_data(data_to_scale):
    scaler = StandardScaler(with_mean=True, with_std=False)
    data_scaled = scaler.fit_transform(data_to_scale)
    scaled_df = pd.DataFrame(data_scaled, index=data_to_scale.index, columns=data_to_scale.columns)
    return scaled_df


def split_data(full_matrix, treatment_labels):
    num_of_samples = full_matrix.shape[0]
    training_size, internal_size = round(num_of_samples * 0.7), round(num_of_samples * 0.3)
    # first separate all 'ambiguous' responders from the training set
    ambigious_samples = [sample_ for sample_ in treatment_labels
                         if treatment_labels[sample_] == 'Ambiguous Responder (101-179 days)']
    ambiguous_set = full_matrix.loc[ambigious_samples, :]
    training_set = full_matrix.drop(ambigious_samples, axis=0)

    # extend internal validation cohort with randomly selected 'good' and 'bad' responders
    num_random_val_samples = internal_size - ambiguous_set.shape[0]
    random_val_samples = training_set.sample(n=num_random_val_samples, random_state=1)
    internal_validation_cohort = pd.concat([ambiguous_set, random_val_samples])
    # remove randomly selected samples from the training set
    training_set = training_set.drop(list(random_val_samples.index.values), axis=0)

    # standard scale TMM input (center data)
    scaled_internal_validation_cohort = standard_scale_data(internal_validation_cohort)
    scaled_training_set = standard_scale_data(training_set)

    # add response labels to matrix
    treatment_labels_internal = {sample_: treatment_labels[sample_]
                                 for sample_ in list(scaled_internal_validation_cohort.index.values)}
    scaled_internal_validation_cohort['response'] = list(scaled_internal_validation_cohort.index.values)
    scaled_internal_validation_cohort['response'] = \
        scaled_internal_validation_cohort['response'].map(treatment_labels_internal)

    treatment_labels_training = {sample_: treatment_labels[sample_]
                                 for sample_ in list(scaled_training_set.index.values)}
    scaled_training_set['response'] = list(scaled_training_set.index.values)
    scaled_training_set['response'] = scaled_training_set['response'].map(treatment_labels_training)

    return scaled_internal_validation_cohort, scaled_training_set


def load_data(path_data, path_metadata):
    metadata = pd.read_csv(path_metadata, delimiter=';')
    metadata = metadata[metadata.matchingRNA == "Yes"]
    treatment_labels = {sample: metadata[metadata['hmfSampleId'] == sample]['Responder'].to_string(index=False).lstrip()
                        for sample in list(metadata['hmfSampleId'])}

    data = pd.read_csv(path_data, delimiter='\t')
    data_transposed = data.T

    return data_transposed, treatment_labels


def run(args):
    input_matrix, treatment_labels = load_data(args.raw_data, args.metadata_path)

    # separate data into training set (70%) and internal validation set (30%)
    #
    internal_validation_cohort, training_set = split_data(input_matrix, treatment_labels)
    internal_validation_cohort.to_csv(os.path.join(args.output_dir, 't_internal_validation_cohort.csv'), header=True,
                                      index=True)
    training_set.to_csv(os.path.join(args.output_dir, 't_training_cohort_scaled.csv'), header=True, index=True)

    spca_model_path, transformed_training_data = run_sparsePCA("training", training_set, int(args.component), args.output_dir)
    transformed_internal_vdata, component_number = transform_with_fitted_sparsePCA(spca_model_path, internal_validation_cohort)
    transformed_internal_vdata.to_csv(
        os.path.join(args.output_dir, "internal_validation_transcriptomics_sparsePCA_{}.csv".format(component_number)),
        header=True, index=True)


args = init()
if args is not 0:
    run(args)