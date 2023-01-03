import csv
import gc
import os
from datetime import datetime
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection as fdr
from IPython.display import display
from attrdict import AttrDict
import matplotlib.pyplot as plt
import numpy as np
import pingouin as pg
from scipy import stats

SRC = "Ubiquitome_AchillesData.xlsx"
NERVOUS_SUFFIX = "CENTRAL_NERVOUS_SYSTEM"
HOMOGENITY_THRESHOLD = 0.5
i = 0

if __name__ == '__main__':
    for NERVOUS_SUFFIX in ['LUNG', 'OVARY', 'ENDOMETRIUM', 'STOMACH', 'CENTRAL_NERVOUS_SYSTEM',
                           'KIDNEY', 'SOFT_TISSUE', 'PANCREAS', 'BONE', 'SKIN', 'BREAST', 'UPPER_AERODIGESTIVE_TRACT',
                           'THYROID', 'LIVER', 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'OESOPHAGUS',
                           'URINARY_TRACT', 'LARGE_INTESTINE', 'AUTONOMIC_GANGLIA']:
        # for NERVOUS_SUFFIX in ['THYROID', 'LIVER', 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'OESOPHAGUS',
        #                      'URINARY_TRACT', 'LARGE_INTESTINE', 'AUTONOMIC_GANGLIA']:
        subjects = ['Total', NERVOUS_SUFFIX, 'Others', 'Delta']  # Gene
        results = list()
        src = pd.read_excel(SRC)
        src = src.reset_index()
        koGenes = AttrDict()
        avgDeltas = AttrDict()
        for _, koGeneData in src.iterrows():
            koGeneKey = koGeneData["Class"] + "-" + koGeneData["Gene"]
            koGeneMeanViability = koGeneData["Mean"]

            koGeneData.pop("index")
            koGeneData.pop("Class")
            koGeneData.pop("Gene")
            koGeneData.pop("Mean")

            koGenes[koGeneKey] = AttrDict({
                "mean": koGeneMeanViability,
                "TtestPValue": None,
                "TT_Type": None,
                "nervous": AttrDict({
                    "values": np.zeros(0),
                    "accumulated_val": 0,
                    "avg": None
                }),
                "notNervous": AttrDict({
                    "values": np.zeros(0),
                    "accumulated_val": 0,
                    "avg": None
                })
            })

            cellLines = koGeneData.keys()
            for cellLine in cellLines:
                is_nervous = "nervous" if NERVOUS_SUFFIX in cellLine.split("_", 1) else "notNervous"
                curr_cell_department = koGenes[koGeneKey][is_nervous]
                curr_cell_department["accumulated_val"] += koGeneData[cellLine]
                curr_cell_department["values"] = np.append(curr_cell_department["values"], koGeneData[cellLine])
            koGenes[koGeneKey]["nervous"]["avg"] = koGenes[koGeneKey]["nervous"]["accumulated_val"] / \
                                                   koGenes[koGeneKey]["nervous"]["values"].size
            koGenes[koGeneKey]["notNervous"]["avg"] = koGenes[koGeneKey]["notNervous"]["accumulated_val"] / \
                                                      koGenes[koGeneKey]["notNervous"]["values"].size
            avgDeltas[koGeneKey] = koGenes[koGeneKey]["nervous"]["avg"] - koGenes[koGeneKey]["notNervous"]["avg"]

            # Assuming Normality assumption (Empiric results)
            # Test Homogeneity assumption:
            is_homogenic = stats.levene(koGenes[koGeneKey]["nervous"]["values"],
                                        koGenes[koGeneKey]["notNervous"]["values"]).pvalue > HOMOGENITY_THRESHOLD
            # Calculate t-test:
            koGenes[koGeneKey]["TT_Type"] = "Student's t-test" if is_homogenic else "Welch's t-test"
            tTestRes = pg.ttest(koGenes[koGeneKey]["nervous"]["values"], koGenes[koGeneKey]["notNervous"]["values"],
                                correction=not is_homogenic)
            koGenes[koGeneKey]["TtestPValue"] = tTestRes["p-val"][0]

            # display(res)

            # ignore irrelevant data:
            if not (koGenes[koGeneKey]["TtestPValue"] < 0.05):
                continue
            # if fdr(np.array([tTestRes["p-val"][0]]), alpha=0.1)[0][0]:
            # print(f"{fdr(np.array([res['p-val'][0]]), alpha=10)[0][0]}")
            # pass
            # print("hi1")
            # continue
            i += 1
            print(i)

            # plot
            scores = [koGenes[koGeneKey]["mean"], koGenes[koGeneKey]["nervous"]["avg"],
                      koGenes[koGeneKey]["notNervous"]["avg"], avgDeltas[koGeneKey]]
            subjects = ['Total', NERVOUS_SUFFIX, 'Others', 'Delta']  # Gene
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.bar(subjects, scores)
            ax.xaxis.set_tick_params(pad=10)
            ax.yaxis.set_tick_params(pad=5)
            ax.grid(b=True, color='gray', linestyle='-.', linewidth=0.75, alpha=0.2)

            # Add Plot Title
            ax.set_title(f'Viability Means - {koGeneKey} KO')

            # Add Text watermark
            fig.text(0.9, 0.15, f"{koGenes[koGeneKey]['TT_Type']} validated", fontsize=12, color='gray', ha='right',
                     va='bottom', alpha=0.7)

            if not os.path.exists(f'Plots/{NERVOUS_SUFFIX}/{koGenes[koGeneKey]["TT_Type"]}'):
                os.makedirs(f'Plots/{NERVOUS_SUFFIX}/{koGenes[koGeneKey]["TT_Type"]}')
            plt.savefig(f'Plots/{NERVOUS_SUFFIX}/{koGenes[koGeneKey]["TT_Type"]}/{koGeneKey} KO.png', dpi=200,
                        format='png',
                        bbox_inches='tight')
            plt.close()
            # scores.append(koGeneKey)
            results.append(
                {
                    "Total": koGenes[koGeneKey]["mean"],
                    NERVOUS_SUFFIX: koGenes[koGeneKey]["nervous"]["avg"],
                    "Others": koGenes[koGeneKey]["notNervous"]["avg"],
                    "Delta": avgDeltas[koGeneKey],
                    "Gene": koGeneKey,
                    "T Test Type": koGenes[koGeneKey]["TT_Type"],
                    "T Test P Value": koGenes[koGeneKey]["TtestPValue"]
                }
            )
        # log results
        subjects.append("Gene")
        subjects.append("T Test Type")
        subjects.append("T Test P Value")
        date = datetime.now().strftime("%Y.%m.%d_%H-%M-%S")
        file_name = f"{date} {NERVOUS_SUFFIX} Viability Score Results.csv"

        with open(file_name, 'w') as results_csv:
            writer = csv.DictWriter(results_csv, subjects)
            writer.writeheader()
            writer.writerows(results)

        csv_data = pd.read_csv(file_name)
        csv_data.to_excel(file_name.replace("csv", "xlsx"), index=None, header=True)
        del results
        del subjects
        gc.collect()

# sorted_deltas = sorted(avgDeltas, key=avgDeltas.get, reverse=True)

# file_name = "{} Averages Sorted Deltas (Big to Small).csv".format(date)
#
# with open(file_name, 'w') as results_csv:
#     writer = csv.writer(results_csv)
#     writer.writerows(zip(avgDeltas.keys(), avgDeltas.values()))
#
# csv_data = pd.read_csv(file_name)
# csv_data.to_excel(file_name.replace("csv", "xlsx"), index=None, header=True)
