import csv
from datetime import datetime

import pandas as pd
from attrdict import AttrDict
import matplotlib.pyplot as plt

SRC = "Ubiquitome_AchillesData.xlsx"
NERVOUS_SUFFIX = "CENTRAL_NERVOUS_SYSTEM"

if __name__ == '__main__':
    src = pd.read_excel(SRC)
    src = src.reset_index()
    koGenes = AttrDict()
    avgDeltas = AttrDict()
    subjects = ['Total', 'Nervous', 'Not Nervous', 'Delta', "Gene"]
    results = [subjects]
    for _, koGeneData in src.iterrows():
        koGeneKey = koGeneData["Class"] + "-" + koGeneData["Gene"]
        koGeneMeanViability = koGeneData["Mean"]

        koGeneData.pop("index")
        koGeneData.pop("Class")
        koGeneData.pop("Gene")
        koGeneData.pop("Mean")

        koGenes[koGeneKey] = AttrDict({
            "mean": koGeneMeanViability,
            "nervous": AttrDict({
                "accumulated_val": 0,
                "amount": 0,
                "avg": None
            }),
            "notNervous": AttrDict({
                "accumulated_val": 0,
                "amount": 0,
                "avg": None
            })
        })

        cellLines = koGeneData.keys()
        for cellLine in cellLines:
            is_nervous = NERVOUS_SUFFIX in cellLine.split("_", 1)
            koGenes[koGeneKey]["nervous" if is_nervous else "notNervous"]["accumulated_val"] += koGeneData[cellLine]
            koGenes[koGeneKey]["nervous" if is_nervous else "notNervous"]["amount"] += 1
            koGenes[koGeneKey]["nervous" if is_nervous else "notNervous"]["avg"] = \
                koGenes[koGeneKey]["nervous" if is_nervous else "notNervous"]["accumulated_val"] / \
                koGenes[koGeneKey]["nervous" if is_nervous else "notNervous"]["amount"]

        avgDeltas[koGeneKey] = koGenes[koGeneKey]["nervous"]["avg"] - koGenes[koGeneKey]["notNervous"]["avg"]
        # plot

        scores = [koGenes[koGeneKey]["mean"], koGenes[koGeneKey]["nervous"]["avg"],
                  koGenes[koGeneKey]["notNervous"]["avg"], avgDeltas[koGeneKey]]
        # fig, ax = plt.subplots(figsize=(7, 5))
        # ax.bar(subjects, scores)
        # ax.xaxis.set_tick_params(pad=10)
        # ax.yaxis.set_tick_params(pad=5)
        # ax.grid(b=True, color='gray', linestyle='-.', linewidth=0.75, alpha=0.2)
        #
        # # Add Plot Title
        # ax.set_title('Viability Means - {} KO'.format(koGeneKey))
        #
        # # Add Text watermark
        # fig.text(0.9, 0.15, "Koren's lab", fontsize=12,
        #          color='gray', ha='right', va='bottom',
        #          alpha=0.7)
        #
        # plt.savefig('Plots/{} KO.png'.format(koGeneKey), dpi=200, format='png', bbox_inches='tight')
        # plt.close()
        scores.append(koGeneKey)
        results.append(scores)

    # log results
    date = datetime.now().strftime("%Y.%m.%d_%H-%M-%S")
    file_name = "{} Neurons Viability Score Results.csv".format(date)

    with open(file_name, 'w') as results_csv:
        writer = csv.writer(results_csv)
        writer.writerows(results)

    csv_data = pd.read_csv(file_name)
    csv_data.to_excel(file_name.replace("csv", "xlsx"), index=None, header=True)

    sorted_deltas = sorted(avgDeltas, key=avgDeltas.get, reverse=True)

    file_name = "{} Averages Sorted Deltas (Big to Small).csv".format(date)

    with open(file_name, 'w') as results_csv:
        writer = csv.writer(results_csv)
        writer.writerows(zip(avgDeltas.keys(), avgDeltas.values()))

    csv_data = pd.read_csv(file_name)
    csv_data.to_excel(file_name.replace("csv", "xlsx"), index=None, header=True)
