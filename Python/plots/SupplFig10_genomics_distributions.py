import seaborn as sns
import matplotlib

## seaborn version: '0.11.2'
## matplotlib version: '3.5.1'
## pandas version: '1.4.2'

matplotlib.style.use('ggplot')
import pandas as pd
from sklearn.preprocessing import StandardScaler

# Values for Genome.TMB", "totalSV", "SV.DUP", "SV.DEL" are available in Suppl. Table 1. - "Mutation Burden" sheet
# For distribution of samples across training and validation sets,
# please refer to "Predictions - HMF" and "Predictions - WCDT" sheets

CPCT_raw_training_genomics = pd.read_csv(".../CPCT_unscaled_training_set_4genomics_w_ccov.csv")
CPCT_raw_training_genomics.set_index("sample", inplace=True)
CPCT_raw_training_genomics = CPCT_raw_training_genomics[["Genome.TMB", "totalSV", "SV.DUP", "SV.DEL"]]
scaler = StandardScaler(with_mean=True, with_std=True)
m_scaled_CPCT_raw_training_genomics = scaler.fit_transform(CPCT_raw_training_genomics.values)
scaled_CPCT_raw_training_genomics = pd.DataFrame(m_scaled_CPCT_raw_training_genomics,
                                                 index=CPCT_raw_training_genomics.index,
                                                 columns=CPCT_raw_training_genomics.columns)
CPCT_raw_training_genomics['cohort'] = ['training' for idx in list(CPCT_raw_training_genomics.index.values)]

CPCT_internal_cohort_genomics = pd.read_csv(".../CPCT_unscaled_internal_validation_cohort_4genomics_w_ccov.csv")
CPCT_internal_cohort_genomics.set_index("sample", inplace=True)
CPCT_internal_cohort_genomics = CPCT_internal_cohort_genomics[["Genome.TMB", "totalSV", "SV.DUP", "SV.DEL"]]
scaler = StandardScaler(with_mean=True, with_std=True)
m_scaled_CPCT_internal_cohort_genomics = scaler.fit_transform(CPCT_internal_cohort_genomics.values)
scaled_CPCT_internal_cohort_genomics = pd.DataFrame(m_scaled_CPCT_internal_cohort_genomics,
                                                 index=CPCT_internal_cohort_genomics.index,
                                                 columns=CPCT_internal_cohort_genomics.columns)
CPCT_internal_cohort_genomics['cohort'] = ['internal' for idx in list(CPCT_internal_cohort_genomics.index.values)]

external_cohort_genomics = pd.read_csv("...")
external_cohort_genomics.set_index("sample", inplace=True)
external_cohort_genomics = external_cohort_genomics[["Genome.TMB", "totalSV", "SV.DUP", "SV.DEL"]]
scaler = StandardScaler(with_mean=True, with_std=True)
m_external_cohort_genomics = scaler.fit_transform(external_cohort_genomics.values)
scaled_external_cohort_genomics = pd.DataFrame(m_external_cohort_genomics,
                                                 index=external_cohort_genomics.index,
                                                 columns=external_cohort_genomics.columns)
external_cohort_genomics['cohort'] = ['external' for idx in list(external_cohort_genomics.index.values)]

# merge raw genomics data from cohorts
all_cohorts = pd.concat([CPCT_raw_training_genomics, CPCT_internal_cohort_genomics, external_cohort_genomics])

scaled_CPCT_internal_cohort_genomics['cohort'] = ['internal' for idx in list(scaled_CPCT_internal_cohort_genomics.index.values)]
scaled_CPCT_raw_training_genomics['cohort'] = ['training' for idx in list(scaled_CPCT_raw_training_genomics.index.values)]
scaled_external_cohort_genomics['cohort'] = ['external' for idx in list(scaled_external_cohort_genomics.index.values)]

# merge scaled genomics data from cohorts
all_cohorts_scaled = pd.concat([scaled_CPCT_raw_training_genomics,
                                scaled_CPCT_internal_cohort_genomics,
                                scaled_external_cohort_genomics])

melted_data = pd.melt(all_cohorts, id_vars=['cohort'], value_vars=["Genome.TMB", "totalSV", "SV.DUP", "SV.DEL"])


with sns.axes_style("white"):
    g = sns.FacetGrid(melted_data, col="variable", hue="cohort", sharex=False, sharey=True, margin_titles=True, height=2.5)
g.map(sns.histplot, "value")
g.add_legend()
g.savefig('.../Suppl_Fig_genomics_distribution_raw.png', dpi=500)


melted_scaled_data = pd.melt(all_cohorts_scaled,
                             id_vars=['cohort'],
                             value_vars=["Genome.TMB", "totalSV", "SV.DUP", "SV.DEL"])

with sns.axes_style("white"):
    m = sns.FacetGrid(melted_scaled_data, col="variable", hue="cohort", sharex=False, sharey=True, margin_titles=True, height=2.5)
m.map(sns.histplot, "value")
m.add_legend()
m.savefig('.../Suppl_Fig_genomics_distribution_scaled.png', dpi=500)