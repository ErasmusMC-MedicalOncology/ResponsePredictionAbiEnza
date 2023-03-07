import pandas as pd
import seaborn as sns
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

matplotlib.style.use('ggplot')

path_suppl_data_table1 = "./Misc/Suppl. Table 1 - OverviewOfData.xls"

coeff_data_df = pd.read_excel(path_suppl_data_table1, sheet_name='Feature importances')
coeff_data_df.rename(columns={"importance": "coefficient"}, inplace=True)
coeff_data_df.replace("Prior_AbiEnza", "prior ARSI", inplace=True)

feature_names = ['Genome.TMB', 'totalSV', 'SV.DUP', 'SV.DEL', 'prior ARSI']

importance_per_fold = []
for i in range(1,80):
    feature_importance = pd.DataFrame(feature_names, columns = ["feature"])
    fold_df = coeff_data_df[coeff_data_df["LOOCV fold"] == i]
    w = np.asarray(list(fold_df["coefficient"]))
    feature_importance["coefficient"] = w
    feature_importance = feature_importance.sort_values(by=["coefficient"], ascending=False)
    importance_per_fold.append(feature_importance)
feature_importance_in_all_folds = pd.concat(importance_per_fold)

sns.set_style("white")
m = sns.stripplot(
    x="coefficient", y="feature", color="#900C3F",
    data=feature_importance_in_all_folds, dodge=True, alpha=0.2
)

sns.barplot(
    data=feature_importance_in_all_folds, x="coefficient", y="feature",
    errcolor='black', errwidth=2,
    errorbar='sd', capsize=.2,
    linewidth=3, edgecolor=".3", facecolor=(0, 0, 0, 0), ax=m
)
m.grid(visible=True, which='major', linestyle='--')

fig = m.get_figure()
plt.tight_layout()
fig.savefig("./Suppl_Fig9_clinicogenomics_feature_importance.png", dpi=500)
