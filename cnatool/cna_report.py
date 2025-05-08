import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

def plot_cna_confusion_matrix_and_summary(
    adata_ground_truth, 
    cna_long, 
    truth_col="simulated_cnvs",
    plot=True
):
    # ---- Parse simulated cnvs ----
    def parse_simulated_cnvs(cnv_string):
        regions = []
        if pd.isna(cnv_string): return regions
        for part in cnv_string.split(','):
            match = re.match(r"(.+?):(\d+)-(\d+)\s*\(CN\s*(-?\d+)\)", part.strip())
            if match:
                chrom, s, e, cnv = match.groups()
                regions.append((str(chrom), int(s), int(e), int(cnv)))
        return regions

    # --- Ground truth table ---
    truth_records = []
    for cell in adata_ground_truth.obs_names:
        for chrom, s, e, cnv in parse_simulated_cnvs(adata_ground_truth.obs.loc[cell, truth_col]):
            if cnv > 2:
                ctype = "gain"
            elif cnv < 2:
                ctype = "loss"
            else:
                ctype = "neutral"
            truth_records.append({"cell": cell, "chrom": str(chrom), "start": int(s), "end": int(e), "true_type": ctype})
    truth_df = pd.DataFrame(truth_records)

    # --- Prediction table ---
    pred_df = cna_long.rename(columns={"chromosome": "chrom", "CNA_call": "pred_type"})
    pred_df = pred_df[["cell", "chrom", "start", "end", "pred_type"]]

    # --- Outer join: include all events seen in either set
    merged = pd.merge(truth_df, pred_df, on=["cell", "chrom", "start", "end"], how="outer")
    merged["true_type"].fillna("neutral", inplace=True)
    merged["pred_type"].fillna("neutral", inplace=True)

    # --- Make confusion matrix ---
    confmat = pd.crosstab(
        merged["true_type"],
        merged["pred_type"],
        rownames=["True"],
        colnames=["Predicted"]
    )

    # --- Calculate summary statistics ---
    classes = ['gain', 'loss', 'neutral']
    summary = {}

    for cls in classes:
        TP = confmat.loc[cls, cls] if cls in confmat.index and cls in confmat.columns else 0
        FP = confmat[cls].sum() - TP if cls in confmat.columns else 0
        FN = confmat.loc[cls].sum() - TP if cls in confmat.index else 0
        precision = TP / (TP + FP) if (TP+FP) > 0 else 0
        recall = TP / (TP + FN) if (TP+FN) > 0 else 0
        f1 = 2*precision*recall/(precision+recall) if precision+recall > 0 else 0
        summary[cls] = {
            "TP": TP, "FP": FP, "FN": FN,
            "precision": precision, "recall": recall, "f1": f1
        }

    # ---- Print summary ----
    print("Class-wise summary:")
    sum_template = "{:<8}{:<8}{:<8}{:<8.2f}{:<8.2f}{:<8.2f}"
    print("{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}".format("CNA", "TP", "FP", "FN", "Prec", "Recall", "F1"))
    for cls in classes:
        s = summary[cls]
        print(sum_template.format(cls, s['TP'], s['FP'], s['FN'], s['precision'], s['recall'], s['f1']))
    print("\nConfusion Matrix:")
    print(confmat)

    if plot:
        plt.figure(figsize=(8, 6))
        sns.heatmap(confmat, annot=True, fmt="d", cmap="Blues")
        plt.title("CNA Confusion Matrix")
        plt.xlabel("Predicted")
        plt.ylabel("True")
        plt.show()

    return confmat, summary

# Usage:
# confmat, summary = plot_cna_confusion_matrix_and_summary(adata_ground_truth, cna_long)