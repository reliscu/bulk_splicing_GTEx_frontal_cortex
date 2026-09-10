import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import rankdata, norm
from statsmodels.stats.multitest import multipletests

# =============================================================================
# 1. SPEARMAN CORRELATION + PERMUTATION P-VALUES + FDR + ANALYTICAL CIs
# =============================================================================

def analytical_ci(r, n, ci=0.95):
    """
    Fisher Z-transform based confidence intervals for correlation.
    Exact for Pearson, approximate for Spearman.

    Parameters
    ----------
    r   : array-like, observed correlations
    n   : int, number of samples
    ci  : float, confidence level (default 0.95)

    Returns
    -------
    lower, upper : arrays of CI bounds
    """
    r = np.asarray(r)
    r_clipped = np.clip(r, -0.9999, 0.9999)
    z = np.arctanh(r_clipped)
    se = 1 / np.sqrt(n - 3)
    z_crit = norm.ppf((1 + ci) / 2)
    lower = np.tanh(z - z_crit * se)
    upper = np.tanh(z + z_crit * se)
    return lower, upper

def spearman_permutation_test(psi, eigengene_df, n_perms=1000, ci=0.95, seed=42):
    """
    Compute Spearman correlations between each PSI event and each cell type,
    with permutation-based p-values, BH-FDR correction, and analytical CIs.
    Also returns permuted correlations for use in permutation-based
    correlation difference test.
 
    Parameters
    ----------
    psi     : pd.DataFrame, shape (n_SEs, n_samples)
    eigengene_df  : pd.DataFrame, shape (n_samples, n_celltypes)
    n_perms         : int, number of permutations
    ci              : float, confidence level for analytical CIs
    seed            : int, random seed
 
    Returns
    -------
    psi_corr_df       : Spearman r values
    psi_pval_df       : permutation p-values
    psi_fdr_df        : BH-FDR corrected p-values (within each cell type)
    psi_ci_lower_df   : analytical CI lower bounds
    psi_ci_upper_df   : analytical CI upper bounds
    perm_corr_results : dict {ct: (n_SEs, n_perms) array} of permuted correlations
    """
    np.random.seed(seed)
 
    n_samples = psi.shape[1]
 
    # --- precompute PSI ranks once ---
    psi_array = psi.values  # (n_SEs, n_samples)
    psi_rank_matrix = np.apply_along_axis(
        lambda x: rankdata(np.where(np.isnan(x), np.nan, x), method='average'),
        axis=1,
        arr=psi_array
    )  # (n_SEs, n_samples)
 
    # center PSI ranks once
    psi_rank_centered = psi_rank_matrix - np.nanmean(psi_rank_matrix, axis=1, keepdims=True)
    psi_rank_norm = np.sqrt(np.nansum(psi_rank_centered ** 2, axis=1, keepdims=True))  # (n_SEs, 1)
 
    # --- precompute all permutation indices ---
    perm_indices = np.array([np.random.permutation(n_samples) for _ in range(n_perms)])
    # shape: (n_perms, n_samples)
 
    corr_results = {}
    pval_results = {}
    perm_corr_results = {}
 
    for ct in eigengene_df.columns:
        ct_vals = eigengene_df[ct].values
        ct_ranks = rankdata(ct_vals, method='average')
 
        # --- observed correlation ---
        ct_centered = ct_ranks - ct_ranks.mean()
        ct_norm = np.sqrt((ct_centered ** 2).sum())
        observed = (psi_rank_centered @ ct_centered) / (psi_rank_norm.squeeze() * ct_norm)
        # shape: (n_SEs,)
 
        # --- all permutations at once ---
        ct_perm_ranks = ct_ranks[perm_indices]          # (n_perms, n_samples)
        ct_perm_centered = ct_perm_ranks - ct_perm_ranks.mean(axis=1, keepdims=True)
        ct_perm_norm = np.sqrt((ct_perm_centered ** 2).sum(axis=1))  # (n_perms,)
 
        # (n_SEs, n_samples) @ (n_samples, n_perms) = (n_SEs, n_perms)
        num = psi_rank_centered @ ct_perm_centered.T
        denom = psi_rank_norm * ct_perm_norm[np.newaxis, :]  # (n_SEs, n_perms)
        perm_corrs = num / denom                             # (n_SEs, n_perms)
 
        # --- two-tailed permutation p-values ---
        pvals = ((np.abs(perm_corrs) >= np.abs(observed[:, np.newaxis])).sum(axis=1) + 1) / (n_perms + 1)
 
        corr_results[ct] = observed
        pval_results[ct] = pvals
        perm_corr_results[ct] = perm_corrs  # store for correlation difference test
 
        print(f"Done: {ct}")
 
    # --- assemble output dataframes ---
    psi_corr_df = pd.DataFrame(corr_results, index=psi.index)
    psi_pval_df = pd.DataFrame(pval_results, index=psi.index)
 
    # --- BH-FDR correction within each cell type ---
    psi_fdr_df = psi_pval_df.copy()
    for ct in eigengene_df.columns:
        pvals = psi_pval_df[ct].values
        mask = ~np.isnan(pvals)
        fdr = np.full(len(pvals), np.nan)
        _, fdr[mask], _, _ = multipletests(pvals[mask], method='fdr_bh')
        psi_fdr_df[ct] = fdr
 
    # --- analytical CIs ---
    n = psi.shape[1]
    psi_ci_lower_df = psi_corr_df.copy()
    psi_ci_upper_df = psi_corr_df.copy()
    for ct in eigengene_df.columns:
        lower, upper = analytical_ci(psi_corr_df[ct].values, n, ci=ci)
        psi_ci_lower_df[ct] = lower
        psi_ci_upper_df[ct] = upper
 
    return (psi_corr_df, psi_pval_df, psi_fdr_df, psi_ci_lower_df, psi_ci_upper_df, perm_corr_results)


# =============================================================================
# 2. PERMUTATION-BASED CORRELATION DIFFERENCE TEST
# =============================================================================

def compare_all_ct_pairs(psi_corr_df, perm_corr_results, ct_cols):
    """
    Permutation-based test for difference between dependent correlations.
    Reuses permuted correlations from spearman_permutation_test

    For each pair (ct1, ct2), the null distribution of r(PSI, ct1) - r(PSI, ct2)
    is computed from the permuted correlations.

    Parameters
    ----------
    psi_corr_df       : pd.DataFrame, observed correlations (no Gene column)
    perm_corr_results : dict {ct: (n_SEs, n_perms) array}
    ct_cols           : list of cell type column names

    Returns
    -------
    dict keyed by (ct1, ct2) tuples, values are DataFrames with columns:
        Gene, r_ct1, r_ct2, r_diff, pval, fdr
    where r_diff = r(PSI, ct1) - r(PSI, ct2), positive means ct1 > ct2
    """
    all_pvals = []
    all_keys = []
    all_masks = []
    results = {}

    for ct1, ct2 in combinations(ct_cols, 2):
        observed_ct1 = psi_corr_df[ct1].values   # (n_SEs,)
        observed_ct2 = psi_corr_df[ct2].values   # (n_SEs,)

        perm_ct1 = perm_corr_results[ct1]         # (n_SEs, n_perms)
        perm_ct2 = perm_corr_results[ct2]         # (n_SEs, n_perms)

        observed_diff = observed_ct1 - observed_ct2           # (n_SEs,)
        perm_diff = perm_ct1 - perm_ct2                       # (n_SEs, n_perms)

        # one-tailed p-value in the observed direction
        pvals = np.where(
            observed_diff >= 0,
            ((perm_diff >= observed_diff[:, np.newaxis]).sum(axis=1) + 1) / (perm_diff.shape[1] + 1),  # p(ct1 > ct2)
            ((perm_diff <= observed_diff[:, np.newaxis]).sum(axis=1) + 1) / (perm_diff.shape[1] + 1) # p(ct1 < ct2)
        )

        mask = ~np.isnan(pvals)
        results[(ct1, ct2)] = {
            'observed_diff': observed_diff,
            'pvals': pvals, 
            'r_ct1': observed_ct1,
            'r_ct2': observed_ct2
        }
        #             'mask': mask,
        
        # we concatenate across all pairs so FDR is corrected globally
        all_pvals.append(pvals[mask])
        all_keys.append((ct1, ct2))
        all_masks.append(mask)

     # --- FDR correct all pairs together, one correction per direction ---
    all_pvals_concat = np.concatenate(all_pvals)
    _, all_fdr_concat, _, _ = multipletests(all_pvals_concat, method='fdr_bh')

    # split FDR back into per-pair results
    idx = 0
    final_results = {}
    for (ct1, ct2), mask in zip(all_keys, all_masks):
        n_valid = mask.sum()

        fdr = np.full(len(mask), np.nan)
        fdr[mask] = all_fdr_concat[idx:idx + n_valid]
        idx += n_valid

        d = results[(ct1, ct2)]
        final_results[(ct1, ct2)] = pd.DataFrame({
            f'r_{ct1}': d['r_ct1'],
            f'r_{ct2}': d['r_ct2'],
            'r_diff': d['observed_diff'],       # positive = ct1 > ct2
            'pval': d['pvals'],  
            'fdr': fdr       
        }, index=psi_corr_df.index)

    return final_results

# =============================================================================
# 3. FIND CELL-TYPE-SPECIFIC SEs
# =============================================================================

def find_specific_SEs_per_ct_basic(psi_fdr_df, psi_corr_df, fdr_thresh=0.05):
    """
    For each cell type, find SEs where FDR < fdr_thresh.

    Returns
    -------
    dict of {ct: DataFrame}
    """
    ct_cols = [c for c in psi_fdr_df.columns if c != 'Gene']
    results = {}
    for target_ct in ct_cols:
        mask = psi_fdr_df[target_ct] < fdr_thresh
        n_total = mask.sum()
        print(f"{target_ct}: {n_total} significant SEs")
        results[target_ct] = psi_corr_df[mask].sort_values(target_ct, ascending=False)
    return results


def find_specific_SEs_per_ct(steiger_results, psi_fdr_df, psi_corr_df, ct_hierarchy, fdr_thresh=0.05, ascending=False):
    """
    For each cell type, find SEs where:
    1. FDR < fdr_thresh for that cell type
    2. Correlation difference is significantly greater (or lesser if ascending)
       than all non-child cell types by permutation-based test (r_diff direction + FDR < fdr_thresh). Children defined in ct_hierarchy are excluded.

    Parameters
    ----------
    steiger_results : dict from compare_all_ct_pairs
    psi_fdr_df      : FDR DataFrame (no Gene column)
    psi_corr_df     : correlation DataFrame (with Gene column)
    fdr_thresh      : float
    ascending       : bool, if True find SEs where target has LOWEST correlation

    Returns
    -------
    dict of {ct: DataFrame}
    """
    
    # ct_hierarchy    : dict {ct: [children to skip]}

    ct_cols = [c for c in psi_fdr_df.columns if c != 'Gene']
    results = {}

    for target_ct in ct_cols:
        children_to_skip = ct_hierarchy.get(target_ct, [])
        # other_cts = [ct for ct in ct_cols if ct != target_ct and ct not in children_to_skip]
        other_cts = [ct for ct in ct_cols if ct != target_ct]
        
        # condition 1: significant in target cell type
        sig_mask = psi_fdr_df[target_ct] < fdr_thresh

        # condition 2: significantly greater (or lesser) than all non-child cell types
        greater_mask = pd.Series(True, index=psi_fdr_df.index)

        for other_ct in other_cts:
            
            if (target_ct, other_ct) in steiger_results:
                df = steiger_results[(target_ct, other_ct)]
                
                if other_ct not in children_to_skip:
                    if ascending:
                        ct_greater = (df['r_diff'] < 0) & (df['fdr'] < fdr_thresh)
                    else:
                        ct_greater = (df['r_diff'] > 0) & (df['fdr'] < fdr_thresh)
                        
                else:
                    # don't require the `r_diff` to be stasticially significant w.r.t. child cell types
                    if ascending:
                        ct_greater = df['r_diff'] < 0
                    else:
                        ct_greater = df['r_diff'] > 0
                        
            elif (other_ct, target_ct) in steiger_results:
                df = steiger_results[(other_ct, target_ct)]
                
                if other_ct not in children_to_skip:
                    if ascending:
                        ct_greater = (df['r_diff'] > 0) & (df['fdr'] < fdr_thresh)
                    else:
                        ct_greater = (df['r_diff'] < 0) & (df['fdr'] < fdr_thresh)
                        
                else:
                    if ascending:
                        ct_greater = df['r_diff'] > 0
                    else:
                        ct_greater = df['r_diff'] < 0
            else:
                print(f"Warning: no result found for ({target_ct}, {other_ct})")
                ct_greater = pd.Series(False, index=psi_fdr_df.index)

            greater_mask = greater_mask & ct_greater

        mask = sig_mask & greater_mask
        n_total = sig_mask.sum()
        n_specific = mask.sum()
        print(f"{target_ct}: {n_specific} specific SEs ({n_total} significant total)")

        results[target_ct] = psi_corr_df[mask].sort_values(target_ct, ascending=ascending)

    return results


# =============================================================================
# 4. COMBINE RESULTS
# =============================================================================

def combine_results(ctype_specific_SEs_basic, ctype_specific_SEs_strict,
                    psi_corr_df, psi_fdr_df, steiger_results, ct_hierarchy):
    """
    For each cell type, combine basic (FDR only) and strict (FDR + permutation
    difference test) results into a single table.

    Steiger columns (r_diff, fdr) are filled only for specific SEs,
    NaN for non-specific SEs.

    Returns
    -------
    dict of {ct: DataFrame}
    """
    ct_cols = [c for c in psi_corr_df.columns if c != 'Gene']
    combined = {}

    for target_ct in ct_cols:
        # start from all FDR-significant SEs
        all_sig = ctype_specific_SEs_basic.get(target_ct, pd.DataFrame()).copy()
        if len(all_sig) == 0:
            continue

        # mark which SEs passed strict criterion in either direction
        strict_true = set(ctype_specific_SEs_strict['True'].get(target_ct, pd.DataFrame()).index)
        strict_false = set(ctype_specific_SEs_strict['False'].get(target_ct, pd.DataFrame()).index)
        all_sig['is_specific'] = all_sig.index.isin(strict_true | strict_false)
        all_sig['specific_direction'] = np.where(
            all_sig.index.isin(strict_false), 'highest',
            np.where(all_sig.index.isin(strict_true), 'lowest', np.nan)
        )

        # add r and fdr for target cell type
        all_sig['r'] = psi_corr_df.loc[all_sig.index, target_ct]
        all_sig['fdr'] = psi_fdr_df.loc[all_sig.index, target_ct]

        # add r_diff and fdr for non-child cell types — only for specific SEs
        children_to_skip = ct_hierarchy.get(target_ct, [])
        other_cts = [ct for ct in ct_cols if ct != target_ct and ct not in children_to_skip]


        for other_ct in other_cts:
            # initialize all as NaN
            all_sig[f'r_diff_{other_ct}'] = np.nan
            all_sig[f'fdr_diff_{other_ct}'] = np.nan

            if (target_ct, other_ct) in steiger_results:
                sdf = steiger_results[(target_ct, other_ct)]
                r_diff_vals = sdf.loc[all_sig.index, 'r_diff'].values
                fdr_vals = sdf.loc[all_sig.index, 'fdr'].values
            elif (other_ct, target_ct) in steiger_results:
                sdf = steiger_results[(other_ct, target_ct)]
                r_diff_vals = r_diff_vals = sdf.loc[all_sig.index,'r_diff'].values * -1  # flip sign
                fdr_vals = fdr_vals = sdf.loc[all_sig.index, 'fdr'].values             # fdr stays the same
            else:
                continue

            all_sig[f'r_diff_{other_ct}'] = r_diff_vals
            all_sig[f'fdr_diff_{other_ct}'] = fdr_vals

        combined[target_ct] = all_sig.sort_values('r', ascending=False)

    return combined
    