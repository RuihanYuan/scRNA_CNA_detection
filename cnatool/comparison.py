import numpy as np
from collections import Counter
from scipy.spatial.distance import pdist
import scanpy as sc
import pandas as pd
from scipy.ndimage import uniform_filter1d
from cnatool.cna_call import add_genome_bins
import scipy

def center_expression(adata_query, adata_ref, region_col="chromosome"):
    """Center gene expression by normal reference mean."""
    means = np.log1p(adata_ref.X).mean(axis=0)
    centered = np.log1p(adata_query.X) - means
    adata_query.X = centered
    return adata_query

def threshold_log_fold_change(adata, max_centered_threshold=0.3):
    """Threshold dynamic range of log-fold-change values."""
    np.clip(adata.X, -max_centered_threshold, max_centered_threshold, out=adata.X)
    return adata

def chromosome_smoothing(adata, region_col="chromosome", window_size=101):
    """Smooth expression values along chromosomes."""
    for chrom in adata.var[region_col].unique():
        gene_indices = np.where(adata.var[region_col] == chrom)[0]
        if len(gene_indices) < window_size:
            continue  # Skip smoothing for short chromosomes
        
        # Apply convolution with uniform filter (pyramidinal weights can be customized if needed)
        adata.X[:, gene_indices] = uniform_filter1d(adata.X[:, gene_indices], size=window_size, axis=1, mode='nearest')
    return adata

def preprocess_adata(adata, adata_ref, region_col="chromosome", max_centered_threshold=3, window_size=101):
    """Preprocess including centering, thresholding, and smoothing."""
    # Center by normal expression
    adata.var = add_genome_bins(adata.var)
    
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
        adata.X = adata.X.astype(np.float32)

    if scipy.sparse.issparse(adata_ref.X):
        adata_ref.X = adata_ref.X.toarray()
        adata_ref.X = adata_ref.X.astype(np.float32)

    adata = center_expression(adata, adata_ref, region_col=region_col)
    
    # Threshold dynamic range
    adata = threshold_log_fold_change(adata, max_centered_threshold=max_centered_threshold)
    
    # Apply chromosome smoothing
    adata = chromosome_smoothing(adata, region_col=region_col, window_size=window_size)
    
    return adata

def cluster_dispersion(adata, cluster_label='cluster', n_pcs=5):
    """Returns dict of mean within-cluster pairwise (PCA-space) distances."""
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata, n_comps=n_pcs)
    disp = {}
    for clust in adata.obs[cluster_label].unique():
        idx = adata.obs[cluster_label] == clust
        X = adata.obsm['X_pca'][idx, :n_pcs]
        if X.shape[0] > 1:
            d = pdist(X)
            disp[clust] = np.mean(d)
        else:
            disp[clust] = 0
    return disp

def cluster_compare(
    adata_ref,
    adata_query,
    cluster_ref_label="cluster",
    cluster_query_label="cluster",
    overlap_threshold=0.5,
    min_cells=10,
    dispersion_ratio_threshold=2.0,
    verbose=True
):
    """
    Compares clusters (by label) from reference and query AnnData objects, flags clusters by 
    group-based properties. No barcode overlap assumed!
    """
    # Ensure same feature set
    genes = list(set(adata_ref.var_names) & set(adata_query.var_names))
    adata_ref = adata_ref[:, genes].copy()
    adata_query = adata_query[:, genes].copy()

    ref_labels = adata_ref.obs[cluster_ref_label].astype(str)
    query_labels = adata_query.obs[cluster_query_label].astype(str)

    # All unique cluster labels in either reference or query
    cluster_names = sorted(list(set(ref_labels.unique()).union(set(query_labels.unique()))))

    dissimilar_clusters = []
    ref_disp = cluster_dispersion(adata_ref, cluster_label=cluster_ref_label)
    query_disp = cluster_dispersion(adata_query, cluster_label=cluster_query_label)

    for cl in cluster_names:
        ref_mask = (ref_labels == cl)
        query_mask = (query_labels == cl)
        ref_n = ref_mask.sum()
        query_n = query_mask.sum()
        if ref_n < min_cells or query_n < min_cells:
            continue  # skip small/unrepresented clusters

        # "Overlap" is now just the ratio of group sizes!
        overlap_rate = min(ref_n, query_n) / max(ref_n, query_n)
        corrected_overlap = overlap_rate  # or, if you want, use your weighted correction scheme

        # Dispersion ratio between query and reference for this cluster label
        this_ref_disp = ref_disp.get(cl, 1e-8)
        this_query_disp = query_disp.get(cl, 1e-8)
        dispersion_ratio = this_query_disp / (this_ref_disp + 1e-8)

        if verbose:
            print(
                f"Cluster/cell_type '{cl}': n_ref={ref_n}, n_query={query_n}, "
                f"overlap_rate={overlap_rate:.2f}, "
                f"dispersion_ratio={dispersion_ratio:.2f}"
            )

        # Flag clusters based on size similarity and/or abnormal dispersion
        if corrected_overlap < overlap_threshold or dispersion_ratio > dispersion_ratio_threshold:
            dissimilar_clusters.append((cl, corrected_overlap, ref_n, query_n, dispersion_ratio))

    # Also return a list of dispersion ratios if needed
    return dissimilar_clusters, [r[-1] for r in dissimilar_clusters]

