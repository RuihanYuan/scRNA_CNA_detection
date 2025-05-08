import scanpy as sc
import pandas as pd
import numpy as np
from collections import Counter

def filter_and_normalize(adata, min_cells_per_gene=10):
    """Filter genes expressed in fewer than specified number of cells and normalize data."""
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    
    # Normalize for sequencing depth (total sum normalization)
    median_sum = np.median(adata.X.sum(axis=1))
    sc.pp.normalize_total(adata, target_sum=median_sum)
    
    return adata

def easy_cluster(
    adata,
    n_pcs=20,
    n_neighbors=10,
    resolution=0.75,
    cluster_label="cluster",
    random_state=0,
    min_cells_per_gene=10
):
    """Standard pipeline with filtering, normalization, and log transformation."""
    
    ad = adata.copy()
    ad.var_names_make_unique()
    ad.obs_names_make_unique()
    
    # Filter and normalize data
    ad = filter_and_normalize(ad, min_cells_per_gene=min_cells_per_gene)
    
    # Log transformation
    sc.pp.log1p(ad)
    # PCA and clustering
    sc.pp.pca(ad, n_comps=n_pcs, random_state=random_state)
    sc.pp.neighbors(ad, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(ad, resolution=resolution, key_added=cluster_label, random_state=random_state)
    
    return ad

def cluster_and_save_reference(adata_ref, out_file, cluster_label="cluster", n_pcs=20, n_neighbors=10, resolution=0.6):
    """Cluster the reference dataset once, save result, and return clustered AnnData."""
    adata_ref = adata_ref.copy()
    adata_ref.var_names_make_unique()
    adata_ref.obs_names_make_unique()
    adata_ref = easy_cluster(adata_ref, n_pcs=n_pcs, n_neighbors=n_neighbors, resolution=resolution, cluster_label=cluster_label)
    adata_ref.write_h5ad(out_file)
    # Optionally, save just the assignments as well
    adata_ref.obs[[cluster_label]].to_csv(out_file.replace(".h5ad", "_clusters.csv"))
    return adata_ref

def save_cluster_labels(adata, cluster_label, out_file):
    adata.obs[[cluster_label]].to_csv(out_file)

def load_cluster_labels(adata, cluster_file, cluster_label):
    # Map saved cluster labels to adata.obs by index/barcode.
    import pandas as pd
    clust = pd.read_csv(cluster_file, index_col=0)
    adata.obs[cluster_label] = clust.loc[adata.obs_names, cluster_label].values