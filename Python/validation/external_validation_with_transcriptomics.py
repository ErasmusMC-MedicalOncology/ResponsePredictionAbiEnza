import os
import argparse
import pandas as pd
import pickle

from precise.pipeline_routine import ConsensusRepresentation

from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

# 'transcriptomics_validation_data_file' format: ENS ids are in the columns,
#                                                and each row is a patient sample
#                                                the ENS ids are sorted in the script so the order does not matter
#

def init():
    parser = argparse.ArgumentParser(prog='transcriptomics_PRECISE.py',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-transcriptomics_ARSI_validation_data_file',
                        help='full file path to transcriptomics + ARSI validation data input')
    parser.add_argument('-transcriptomics_training_data_file',
                        help='full file path to "CPCT_cohort_NEW_transcriptomics_training_set.csv"')
    parser.add_argument('-genomics_ARSI_validation_data_file',
                        help='full file path to the genomics + ARSI validation data input')
    parser.add_argument('-genomics_ARSI_training_data_file',
                        help='full file path to "CPCT_cohort_genomics_ARSI_data.csv"')
    parser.add_argument('-genomics_model_file', help='full path of genomics.pkl')
    parser.add_argument('-output_folder',
                        help='output folder name where prediction output and intermedier files are also saved')

    args = parser.parse_args()
    return args


def load_model(model_file_path):
    with open(model_file_path, 'rb') as file:
        model_ = pickle.load(file)
    return model_


def train_linSVC_classifier(X_training, y_training):
    ### train calibrated linSVC classifier
    l_svc = LinearSVC(dual=False, random_state=1, max_iter=10000)
    ccf = CalibratedClassifierCV(l_svc)
    ccf.fit(X_training, y_training)

    return ccf


def run(args):
    ### genomics + ARSI features
    genomic_features = ['Genome.TMB', 'totalSV', 'SV.DUP', 'SV.DEL']
    clinvar_features = ['Prior_AbiEnza']

    # load data
    training_set_with_response = pd.read_csv(args.transcriptomics_training_data_file, index_col=[0])
    validation_set_unsorted = pd.read_csv(args.transcriptomics_ARSI_validation_data_file, index_col=[0])
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

    # subset the dataframe with the Prior_AbiEnza variable from the transcriptomics + Prior_AbiEnza data
    X_valid_clinvar_transcriptomics = validation_set_unsorted[clinvar_features]
    # Remove NAs. These are also the same as with full clincal FU, so we do not lose any samples in the final anlaysis.
    X_valid_clinvar_transcriptomics.dropna(inplace=True)

    ### standard scale data
    X_training_data = training_set.values
    X_validation_data = validation_set.values

    ### consensus computation
    consensus = ConsensusRepresentation(source_data=X_training_data,
                                        target_data=X_validation_data,
                                        n_factors=40,
                                        n_pv=40,
                                        dim_reduction='ica',
                                        n_representations=40,
                                        use_data=False,
                                        mean_center=True,
                                        std_unit=True)

    consensus.fit(X_training_data)

    # save consensus model
    pickle.dump(consensus,
                open(os.path.join(args.output_folder, "consensus_model_40_ICs.pkl"), 'wb'))

    X_training_transformed_ = consensus.transform(X_training_data)
    X_validation_transformed_ = consensus.transform(X_validation_data)

    X_training_transformed = pd.DataFrame(X_training_transformed_,
                columns=["C{}".format(c) for c in range(40)],
                index=training_set.index)
    X_validation_transformed = pd.DataFrame(X_validation_transformed_,
                columns=["C{}".format(c) for c in range(40)],
                index=validation_set.index)

    # load validation genomics + ARSI data
    valid_genomics_arsi_df = pd.read_csv(args.genomics_ARSI_validation_data_file, index_col=[0])

    X_valid_gen_ = valid_genomics_arsi_df[genomic_features]
    # prescale genomics features (only validation)
    scaler_pregen = StandardScaler(with_mean=True, with_std=True)
    X_valid_gen = pd.DataFrame(
        scaler_pregen.fit_transform(X_valid_gen_),
        columns=X_valid_gen_.columns,
        index=X_valid_gen_.index)

    # Prior_AbiEnza: X_valid_clinvar_transcriptomics, X_valid_clinvar_genomics
    X_valid_clinvar_genomics = valid_genomics_arsi_df[clinvar_features]
    X_valid_gen_clinvar = X_valid_gen.join(X_valid_clinvar_genomics)

    X_valid_gen_clinvar_transcriptomics = X_valid_gen_clinvar.join(X_validation_transformed)
    X_valid_clinvar_transcriptomics = X_valid_clinvar_transcriptomics.join(X_validation_transformed)

    ### center validation data
    scaler_gen_ARSI = StandardScaler(with_mean=True, with_std=False)
    ##
    X_valid_gen_clinvar_transcriptomics_scaled = pd.DataFrame(scaler_gen_ARSI.fit_transform(X_valid_gen_clinvar_transcriptomics),
                                columns=X_valid_gen_clinvar_transcriptomics.columns,
                                index=X_valid_gen_clinvar_transcriptomics.index)
    ## the 'X_valid_gen_clinvar_transcriptomics_scaled' still has the 'Prior_AbiEnza' from the genomics data file
    ## BUT the WGS-only sample with 'NaN' transcriptomics features will be removed here:
    X_valid_gen_clinvar_transcriptomics_scaled.dropna(inplace=True)
    ## get the list of sample IDs used for the validation of the transcriptomics + genomics + ARSI model
    validation_samples_transcriptomics_gen_ARSI = list(X_valid_gen_clinvar_transcriptomics_scaled.index.values)

    scaler_ARSI = StandardScaler(with_mean=True, with_std=False)
    ## the 'X_valid_clinvar_transcriptomics_scaled' now should have the 'Prior_AbiEnza' parsed from the transcriptomics data file
    X_valid_clinvar_transcriptomics_scaled = pd.DataFrame(
        scaler_ARSI.fit_transform(X_valid_clinvar_transcriptomics),
        columns=X_valid_clinvar_transcriptomics.columns,
        index=X_valid_clinvar_transcriptomics.index)
    ## get the list of sample IDs used for the validation of the transcriptomics + ARSI model
    validation_samples_transcriptomics_ARSI = list(X_valid_clinvar_transcriptomics_scaled.index.values)


    # load training genomics + ARSI data
    train_genomics_arsi_df = pd.read_csv(args.genomics_ARSI_training_data_file, index_col=[0])

    X_train_gen = train_genomics_arsi_df[genomic_features]
    X_train_clinvar = train_genomics_arsi_df[clinvar_features]
    X_train_gen_clinvar = X_train_gen.join(X_train_clinvar)

    X_train_gen_clinvar_transcriptomics = X_train_gen_clinvar.join(X_training_transformed)
    X_train_clinvar_transcriptomics = X_train_clinvar.join(X_training_transformed)

    ### center training data
    scaler_gen_ARSI = StandardScaler(with_mean=True, with_std=False)
    X_train_gen_clinvar_transcriptomics_scaled = pd.DataFrame(
        scaler_gen_ARSI.fit_transform(X_train_gen_clinvar_transcriptomics),
        columns=X_train_gen_clinvar_transcriptomics.columns,
        index=X_train_gen_clinvar_transcriptomics.index)

    scaler_ARSI = StandardScaler(with_mean=True, with_std=False)
    X_train_clinvar_transcriptomics_scaled = pd.DataFrame(
        scaler_ARSI.fit_transform(X_train_clinvar_transcriptomics),
        columns=X_train_clinvar_transcriptomics.columns,
        index=X_train_clinvar_transcriptomics.index)

    #########################################
    ### train transcriptomics-only classifier
    trained_model = train_linSVC_classifier(X_training_transformed, y_training)
    # save model
    pickle.dump(trained_model,
                open(os.path.join(args.output_folder, "trained_40_ICs_PRECISE_based_model.pkl"), 'wb'))

    ### predict validation cohort labels
    test_prob_tr_only = trained_model.predict_proba(X_validation_transformed)
    pred_label = trained_model.predict(X_validation_transformed)
    test_p = test_prob_tr_only

    prediction_output_file = os.path.join(args.output_folder,  "predictions_consensus_model_40_ICs.csv")
    prediction_output = open(prediction_output_file, 'w')
    prediction_output.write("sampleID,predicted_label,class_0_probability,class_1_probability\n")

    for vx, val_sample in enumerate(validation_samples):
        sample_test_p = test_p[vx]
        prediction_output.write(
            "{},{},{},{}\n".format(val_sample, str(pred_label[vx]), str(sample_test_p[0]), str(sample_test_p[1])))
    prediction_output.close()

    #########################################
    ### train transcriptomics + ARSI classifier
    trained_model = train_linSVC_classifier(X_train_clinvar_transcriptomics_scaled, y_training)
    # save model
    pickle.dump(trained_model,
                open(os.path.join(args.output_folder, "trained_40_ICs_ARSI_PRECISE_based_model.pkl"), 'wb'))

    ### predict validation cohort labels
    test_prob = trained_model.predict_proba(X_valid_clinvar_transcriptomics_scaled)
    pred_label = trained_model.predict(X_valid_clinvar_transcriptomics_scaled)
    test_p = test_prob

    prediction_output_file = os.path.join(args.output_folder,  "predictions_consensus_model_40_ICs_ARSI.csv")
    prediction_output = open(prediction_output_file, 'w')
    prediction_output.write("sampleID,predicted_label,class_0_probability,class_1_probability\n")

    for vx, val_sample in enumerate(validation_samples_transcriptomics_ARSI):
        sample_test_p = test_p[vx]
        prediction_output.write(
            "{},{},{},{}\n".format(val_sample, str(pred_label[vx]), str(sample_test_p[0]), str(sample_test_p[1])))
    prediction_output.close()

    #########################################
    ### train transcriptomics + genomics + ARSI classifier
    trained_model = train_linSVC_classifier(X_train_gen_clinvar_transcriptomics_scaled, y_training)
    # save model
    pickle.dump(trained_model,
                open(os.path.join(args.output_folder, "trained_40_ICs_genomics_ARSI_PRECISE_based_model.pkl"), 'wb'))

    ### predict validation cohort labels
    test_prob = trained_model.predict_proba(X_valid_gen_clinvar_transcriptomics_scaled)
    pred_label = trained_model.predict(X_valid_gen_clinvar_transcriptomics_scaled)
    test_p = test_prob

    prediction_output_file = os.path.join(args.output_folder,  "predictions_consensus_model_40_ICs_ARSI_genomics.csv")
    prediction_output = open(prediction_output_file, 'w')
    prediction_output.write("sampleID,predicted_label,class_0_probability,class_1_probability\n")

    for vx, val_sample in enumerate(validation_samples_transcriptomics_gen_ARSI):
        sample_test_p = test_p[vx]
        prediction_output.write(
            "{},{},{},{}\n".format(val_sample, str(pred_label[vx]), str(sample_test_p[0]), str(sample_test_p[1])))
    prediction_output.close()

    #########################################
    ### Averaging ensemble
    ### 1) run predictions with genomics model
    ### 2) calculate averaging ensemble result - genomics + transcriptomics (40 ICs)

    # for the genomics model predictions, the same sample set is used as for the transcriptomics + ARSI + genomics model
    # which is why the 'X_valid_gen_clinvar_transcriptomics' dataframe is used here
    genomics_model = load_model(args.genomics_model_file)
    ### Make predictions
    X_valid_gen_clinvar_transcriptomics.dropna(inplace=True)
    X_genomics_only = X_valid_gen_clinvar_transcriptomics[genomic_features]
    gen_prob = genomics_model.predict_proba(X_genomics_only)

    # take predicted probabilities only for samples that are both in the genomics and transcriptomics validation
    validation_samples_gen = list(X_genomics_only.index.values)
    validation_samples_tr = list(X_validation_transformed.index.values) # transcriptomics-only validation samples
    common_val_samples = list(set(validation_samples_gen) & set(validation_samples_tr))

    ### save predictions
    prediction_output_file_avg_ens = os.path.join(args.output_folder,  "predictions_averaging_ensemble.csv")
    prediction_output_avg_ens = open(prediction_output_file_avg_ens, 'w')
    prediction_output_avg_ens.write("sampleID,predicted_label,class_0_probability,class_1_probability\n")

    for val_sample in common_val_samples:
        # get sample index
        gen_idx = validation_samples_gen.index(val_sample)
        tr_idx = validation_samples_tr.index(val_sample)

        p_gen = gen_prob[gen_idx] #predicted p from genomics-only
        p_tr = test_prob_tr_only[tr_idx] #predicted p from transcriptomics-only model

        class_1_p = (float(p_gen[1]) + float(p_tr[1])) / 2
        class_0_p = (float(p_gen[0]) + float(p_tr[0])) / 2
        ae_pred_label = 1 if class_1_p > class_0_p else 0

        prediction_output_avg_ens.write(
            "{},{},{},{}\n".format(val_sample, str(ae_pred_label), str(class_0_p), str(class_1_p)))
    prediction_output_avg_ens.close()




args = init()
if args is not 0:
    run(args)