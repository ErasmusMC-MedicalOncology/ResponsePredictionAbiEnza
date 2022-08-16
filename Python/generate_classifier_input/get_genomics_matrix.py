import os
import argparse
import pandas as pd
from sklearn.preprocessing import StandardScaler

def init():
    parser = argparse.ArgumentParser(prog='train_genomics_arsi.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-metadata_path', help='metadata file path')
    parser.add_argument('-raw_data', help='raw data path')
    parser.add_argument('-output_dir', help='output directory')
    args = parser.parse_args()

    if not args.metadata_path or not args.raw_data or not args.output_dir:
        parser.print_help()
        return 0
    return args


def standard_scale_data(data_to_scale):
    scaler = StandardScaler(with_mean=True, with_std=True)
    data_scaled = scaler.fit_transform(data_to_scale)
    scaled_df = pd.DataFrame(data_scaled, index=data_to_scale.index, columns=data_to_scale.columns)
    return scaled_df


def include_clinical_covariates(gtraining_set, metadata_path):
    gtraining_set_samples = list(gtraining_set.index.values)
    metadata = pd.read_csv(metadata_path, delimiter=';')
    metadata = metadata[metadata.hmfSampleId.isin(gtraining_set_samples)]
    metadata['Prior_AbiEnza'] = metadata[['Prior_Abiraterone', 'Prior_Enzalutamide']].max(axis=1)
    metadata['Prior_Chemo'] = metadata[['Prior_Otherchemotherapy', 'Prior_Docetaxel', 'Prior_Cabazitaxel']].max(axis=1)


    ntlines_labels = {
        sample: metadata[metadata['hmfSampleId'] == sample]['Number_prior_treatment_lines'].to_string(
            index=False).lstrip().replace(",", ".")
        for sample in list(gtraining_set.index.values)}
    chemo_labels = {
        sample: metadata[metadata['hmfSampleId'] == sample]['Prior_Chemo'].to_string(index=False).lstrip()
        for sample in list(gtraining_set.index.values)}
    arsi_labels = {
        sample: metadata[metadata['hmfSampleId'] == sample]['Prior_AbiEnza'].to_string(index=False).lstrip()
        for sample in list(gtraining_set.index.values)}
    gtraining_set['Number_prior_treatment_lines'] = list(gtraining_set.index.values)
    gtraining_set['Number_prior_treatment_lines'] = gtraining_set['Number_prior_treatment_lines'].map(ntlines_labels)

    gtraining_set['Prior_Chemo'] = list(gtraining_set.index.values)
    gtraining_set['Prior_Chemo'] = gtraining_set['Prior_Chemo'].map(chemo_labels)

    gtraining_set['Prior_AbiEnza'] = list(gtraining_set.index.values)
    gtraining_set['Prior_AbiEnza'] = gtraining_set['Prior_AbiEnza'].map(arsi_labels)

    gtraining_set = gtraining_set.replace(to_replace=['Yes', 'No'], value=[1, 0])


    return gtraining_set


def load_genomics_data(path_gdata, path_metadata, training_samples):
    metadata = pd.read_csv(path_metadata, delimiter=';')
    metadata = metadata[metadata.hmfSampleId.isin(training_samples)]
    treatment_labels_internal = {
        sample: metadata[metadata['hmfSampleId'] == sample]['Responder'].to_string(index=False).lstrip()
        for sample in list(metadata['hmfSampleId']) if sample in training_samples}

    gdata = pd.read_excel(path_gdata, 'Mutation Burden')
    gdata = gdata[['Genome.TMB', 'totalSV', 'SV.DUP', 'SV.DEL', 'sample']]
    gdata = gdata[gdata['sample'].isin(training_samples)]
    gdata = gdata.set_index('sample')

    scaled_genomics_data = standard_scale_data(gdata)
    scaled_genomics_data['response'] = list(scaled_genomics_data.index.values)
    scaled_genomics_data['response'] = scaled_genomics_data['response'].map(treatment_labels_internal)

    return scaled_genomics_data


def run(args):
    # get separated data (training set (70%) and internal validation set (30%))
    # based on the already analysed transcriptomics dataset
    #
    internal_validation_cohort = pd.read_csv(
        args.internal_cohort,
        index_col=[0])
    training_set = pd.read_csv(
        args.training_cohort,
        index_col=[0])
    training_set_samples = list(training_set.index.values)
    validation_set_samples = list(internal_validation_cohort.index.values)

    # genomic features are scaled and centered
    #
    training_cohort_genomics_matrix = load_genomics_data(args.raw_data, args.metadata_path, training_set_samples)
    internal_validation_cohort_genomics_matrix = load_genomics_data(args.raw_data, args.metadata_path, validation_set_samples)

    training_cohort_genomics_matrix.to_csv(os.path.join(args.output_dir, "training_set_4genomics.csv"), header=True, index=True)
    internal_validation_cohort_genomics_matrix.to_csv(os.path.join(args.output_dir, "internal_validation_set_4genomics.csv"), header=True, index=True)

    # add clinical covariates to the matrices
    #
    training_cohort_genomics_matrix_with_ccovs = include_clinical_covariates(training_cohort_genomics_matrix, args.metadata_path)
    internal_validation_cohort_genomics_matrix_with_ccovs = include_clinical_covariates(internal_validation_cohort_genomics_matrix, args.metadata_path)
    training_cohort_genomics_matrix_with_ccovs.to_csv(
        os.path.join(args.output_dir, 'training_set_4genomics_w_ccov.csv'),
        header=True, index=True)
    internal_validation_cohort_genomics_matrix_with_ccovs.to_csv(
        os.path.join(args.output_dir, 'internal_validation_set_4genomics_w_ccov.csv'),
        header=True, index=True)

args = init()
if args is not 0:
    run(args)
