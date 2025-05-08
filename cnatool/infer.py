import os
import numpy as np
import pandas as pd
import scanpy as sc
from cnatool.clustering import easy_cluster
from cnatool.comparison import cluster_compare, preprocess_adata
from cnatool.cna_call import call_cna_by_region, add_genome_bins
import scipy.sparse


DEFAULT_REFERENCEDIR = r"D:\bme graduate\EN.580.447\Project Package 2\cnatool\reference"
DEFAULT_CLUSTERED = os.path.join(DEFAULT_REFERENCEDIR, "reference_clustered.h5ad")
DEFAULT_UNCLUSTERED = os.path.join(DEFAULT_REFERENCEDIR, "reference.h5ad")

def ensure_clustered_reference(clustered_path, unclustered_path, cluster_label="cluster"):
    if os.path.exists(clustered_path):
        return clustered_path
    elif os.path.exists(unclustered_path):
        print(f"{clustered_path} not found, generating from {unclustered_path}...")
        from .clustering import cluster_and_save_reference
        adata_ref = sc.read_h5ad(unclustered_path)
        adata_ref.var_names_make_unique()
        adata_ref.obs_names_make_unique()
        cluster_and_save_reference(adata_ref, out_file=clustered_path, cluster_label=cluster_label)
        return clustered_path
    else:
        raise FileNotFoundError(
            f"Neither clustered ({clustered_path}) nor unclustered ({unclustered_path}) reference found."
        )

def build_cna_long_format(cna_matrix, adata_var, cluster_id, cna_types=('gain', 'loss')):
    rows = []
    bins = [c for c in cna_matrix.columns if c in adata_var['region'].values]
    for cell in cna_matrix.index:
        for reg in bins:
            call = cna_matrix.at[cell, reg]
            if call in cna_types:
                meta = adata_var.loc[adata_var['region'] == reg]
                if not meta.empty:
                    chromosome = meta['chromosome'].iloc[0]
                    start = meta['start'].iloc[0]
                    end = meta['end'].iloc[0]
                else:
                    chromosome = start = end = None
                rows.append({
                    'cell': cell,
                    'chromosome': chromosome,
                    'start': start,
                    'end': end,
                    'region': reg,
                    'CNA_call': call,
                    'cluster': cluster_id
                })
    return pd.DataFrame(rows)

def center_cells(adata):
    medians = np.median(adata.X, axis=1)
    adata.X = adata.X - medians[:, np.newaxis]
    return adata

def infer_cna(
    adata,
    adata_ref=None,
    ref_clustered_file=None,
    cluster_label="cluster",
    overlap_threshold=0.5,
    min_cells=10,
    dispersion_ratio_threshold=2.0,
    min_genes=5,
    z_threshold=1.0,
    region_col="region",             # make region_col explicit
    max_centered_threshold=3,
    window_size=101,
    verbose=True
):
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Load or prepare reference
    if adata_ref is None and ref_clustered_file is None:
        default_ref_path = ensure_clustered_reference(DEFAULT_CLUSTERED, DEFAULT_UNCLUSTERED, cluster_label)
        if verbose:
            print(f"No adata_ref/ref_clustered_file supplied, using default reference: {default_ref_path}")
        from .io import load_clustered_reference
        adata_ref_clust = load_clustered_reference(default_ref_path, cluster_label)
    elif adata_ref is not None:
        if cluster_label in adata_ref.obs.columns:
            adata_ref_clust = adata_ref
        else:
            adata_ref_clust = easy_cluster(adata_ref, cluster_label=cluster_label)
    elif ref_clustered_file is not None:
        from .io import load_clustered_reference
        adata_ref_clust = load_clustered_reference(ref_clustered_file, cluster_label)
    else:
        raise ValueError("Provide adata_ref or ref_clustered_file.")

    # Cluster query if needed
    adata_query_clust = easy_cluster(adata, cluster_label=cluster_label)

    # Center query cells
    adata_query_clust = center_cells(adata_query_clust)

    # Adjust relative to normal reference cells
    common_genes = adata_query_clust.var_names.intersection(adata_ref_clust.var_names)
    adata_query_clust = adata_query_clust[:, common_genes].copy()
    adata_ref_clust = adata_ref_clust[:, common_genes].copy()
    ref_means = adata_ref_clust.X.mean(axis=0)
    adata_query_clust.X = adata_query_clust.X - ref_means

    # Revert log transformation
    adata_query_clust.X = np.expm1(adata_query_clust.X)

    


    # 1. Compare clusters and flag abnormal groups
    disclusters, reasons = cluster_compare(
        adata_ref_clust,
        adata_query_clust,
        cluster_label,
        cluster_label,
        overlap_threshold=overlap_threshold,
        min_cells=min_cells,
        dispersion_ratio_threshold=dispersion_ratio_threshold,
        verbose=verbose
    )

    if not disclusters:
        print("No abnormal clusters found at given thresholds.")
        return disclusters, pd.DataFrame(), pd.DataFrame()

    # ----- MAIN LOOP with cluster-based preprocessing -----
    all_summaries = []
    all_long_assignments = []
    for cl, reason in zip(disclusters, reasons):
        cluster_id = cl[0]
        if verbose:
            print(f"Flagged cluster {cluster_id}: {reason} -- Calling CNAs for this cluster...")

        # Select cells in this abnormal query cluster
        query_mask = (adata_query_clust.obs[cluster_label] == cluster_id)
        adata_sub_query = adata_query_clust[query_mask].copy()
        if scipy.sparse.issparse(adata_sub_query.X):
            adata_sub_query.X = adata_sub_query.X.toarray()
        adata_sub_query.X = adata_sub_query.X.astype(np.float32)

        # Select corresponding cells in the reference cluster
        ref_mask = (adata_ref_clust.obs[cluster_label] == cluster_id)
        adata_sub_ref = adata_ref_clust[ref_mask].copy()
        if scipy.sparse.issparse(adata_sub_ref.X):
            adata_sub_ref.X = adata_sub_ref.X.toarray()
        adata_sub_ref.X = adata_sub_ref.X.astype(np.float32)

        # ---- PREPROCESS per cluster (the KEY step) ----
        proc_query = preprocess_adata(
            adata_sub_query.copy(),
            adata_sub_ref.copy(),
            region_col=region_col,
            max_centered_threshold=max_centered_threshold,
            window_size=window_size
        )
        proc_ref = preprocess_adata(
            adata_sub_ref.copy(),
            adata_sub_ref.copy(),
            region_col=region_col,
            max_centered_threshold=max_centered_threshold,
            window_size=window_size
        )
        # -----------------------------------------------

        # --- Per-cell-by-bin calling (on preprocessed) ---
        cna_calls = call_cna_by_region(
            proc_query,
            reference_cells=proc_ref,   # This can be AnnData or boolean mask, as your function expects
            min_genes=min_genes,
            z_threshold=z_threshold
        )

        # Build long-format CNA assignment table
        long_df = build_cna_long_format(
            cna_calls.drop(columns=['cluster'], errors='ignore'),
            proc_query.var,
            cluster_id
        )
        all_long_assignments.append(long_df)

        # Prepare per-cluster summary DataFrame
        cna_calls['cluster'] = cluster_id
        calls_long = (
            cna_calls
            .reset_index()
            .melt(id_vars=['index', 'cluster'], var_name='region', value_name='call')
            .rename(columns={'index': 'cell'})
        )

        summary = (
            calls_long
            .groupby(['cluster', 'region', 'call'])
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )

        # Join with region info
        gene_annot = proc_query.var.reset_index().rename(columns={'index': 'region'})
        summary = summary.merge(
            gene_annot[['region', 'chromosome', 'start', 'end', 'strand']],
            on='region',
            how='left'
        )
        col_order = ['cluster', 'region', 'chromosome', 'start', 'end', 'strand'] + [
            c for c in summary.columns if c in {'gain', 'loss', 'neutral'}
        ]
        summary = summary[[c for c in col_order if c in summary.columns]]
        all_summaries.append(summary)
        if verbose:
            print(summary.head())

    # Concatenate all abnormal clusters' results
    final_summary = pd.concat(all_summaries, ignore_index=True) if all_summaries else pd.DataFrame()
    full_cna_long = pd.concat(all_long_assignments, ignore_index=True) if all_long_assignments else pd.DataFrame()

    # Add to AnnData.uns for downstream use
    adata_query_clust.uns['cna_summary'] = final_summary
    adata_query_clust.uns['cna_cell_long'] = full_cna_long

    return disclusters, final_summary, full_cna_long