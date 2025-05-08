import numpy as np
import pandas as pd

def add_genome_bins(var, bin_size=1000000):
    """Assigns each gene to a bin based on chromosome and start coordinate."""
    var = var.copy()
    needed = {'chromosome', 'start', 'end'}
    assert needed.issubset(var.columns), f"Missing columns in .var: {needed - set(var.columns)}"
    nonna = var.dropna(subset=['chromosome', 'start'])
    nonna['start'] = nonna['start'].astype(int)
    # Calculate bin start
    nonna['bin_start'] = (nonna['start'] // bin_size) * bin_size
    nonna['bin_end']   = nonna['bin_start'] + bin_size - 1
    nonna['region'] = (nonna['chromosome'].astype(str)
                       + ':' + nonna['bin_start'].astype(str)
                       + '-' + nonna['bin_end'].astype(str))
    # Merge labels back into .var (fill missing as needed)
    var['region'] = nonna['region']
    var['bin_start'] = nonna['bin_start']
    var['bin_end'] = nonna['bin_end']
    return var


def call_cna_by_region(adata, reference_cells=None, min_genes=10, z_threshold=1.5):
    """
    Call gain/loss by region (bin) for each cell by comparing mean expr to reference.
    Adds DataFrame to adata.uns['cna_calls'].
    Regions are determined by adata.var['region']!
    """
    if reference_cells is None:
        reference_cells = np.ones(adata.n_obs, dtype=bool)

    gene_regions = adata.var['region']
    bins = gene_regions.dropna().unique()

    X = adata.raw.X if adata.raw is not None else adata.X
    if not isinstance(X, np.ndarray): X = X.toarray()
    df_expr = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)

    results = pd.DataFrame(index=adata.obs_names, columns=bins)
    for reg in bins:
        genes = gene_regions[gene_regions == reg].index
        if len(genes) < min_genes: continue
        region_expr = df_expr[genes].mean(axis=1)
        if isinstance(reference_cells, np.ndarray) and reference_cells.dtype == bool:
            ref_mean = region_expr[reference_cells].mean()
            ref_std = region_expr[reference_cells].std()
        elif isinstance(reference_cells, (pd.Series, list, np.ndarray)):
            mask = region_expr.index.isin(reference_cells)
            ref_mean = region_expr[mask].mean()
            ref_std = region_expr[mask].std()
        else:
            # fallback: use all cells
            ref_mean = region_expr.mean()
            ref_std = region_expr.std()
        zscore = (region_expr - ref_mean)/(ref_std + 1e-8)
        calls = pd.Series('neutral', index=region_expr.index)
        calls[zscore > z_threshold] = 'gain'
        calls[zscore < -z_threshold] = 'loss'
        results[reg] = calls
    adata.uns['cna_calls'] = results
    return results