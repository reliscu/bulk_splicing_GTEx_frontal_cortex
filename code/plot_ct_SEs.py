import re
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

def _safe(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name))

def plot_corr_distros(
    psi_corr_df,
    psi_fdr_df,
    data_source,
    outdir,
    fdr_thresh=0.05
):
    
    ct_cols = [c for c in psi_fdr_df.columns if c != 'Gene']
    os.makedirs(outdir, exist_ok=True)
    
    n_cols = 2
    n_rows = int(np.ceil(len(ct_cols) / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 2.5 * n_rows), dpi=200, sharex=True) 

    axes_flat = axes.flatten()

    for i, ct in enumerate(ct_cols):
        ax = axes_flat[i]
        observed = psi_corr_df[ct].values
        sig_mask = (psi_fdr_df[ct] < fdr_thresh).values

        ax.tick_params(labelbottom=True)
        ax.hist(observed, bins=100, color='steelblue', alpha=0.6, density=False, label='Observed (all)')
        ax.hist(observed[sig_mask], bins=100, color='salmon', alpha=0.8, density=False, 
                label=f'Significant (FDR < {fdr_thresh})')
        ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
        ax.set_title(ct)
        ax.set_ylabel('Density')
        ax.set_xlabel('Spearman r')
        ax.legend(fontsize=7)

    # hide any unused axes
    for j in range(len(ct_cols), len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.tight_layout()
    plt.savefig(f"{outdir}/{data_source}_corr_distributions.pdf", dpi=200)
    plt.close()

def plot_ct_SE_barplots(
    ct_SEs,
    psi_corr_df,
    psi_fdr_df,
    expr_corr_df,
    data_source,
    outdir,
    psi_ci_lower_df,   
    psi_ci_upper_df,
    ct_cols=None,
    fdr_thresh=0.05,
    specific=True,
    direction="highest",
    top_n=25
):
    os.makedirs(outdir, exist_ok=True)
    
    if not ct_cols:
        ct_cols = [c for c in psi_corr_df.columns if c != 'Gene']

    for target_ct, se_df in ct_SEs.items():
        
        if specific:
            se_df_subset = se_df[se_df['specific_direction'] == direction]

            if len(se_df_subset) == 0:
                print(f"No {target_ct} specific SEs.")
                continue

            print(f"Plotting {len(se_df_subset)} {target_ct} specific SEs.")
                
        else:
            other_cts = [ct for ct in ct_cols if ct != target_ct]

            if direction == "highest":
                se_df = se_df[se_df['r'] > 0]
                mask = se_df[target_ct] > se_df[other_cts].max(axis=1)
                
            else:
                se_df = se_df[se_df['r'] < 0].sort_values(target_ct)
                mask = se_df[target_ct] < se_df[other_cts].min(axis=1)
            
            n = min(sum(mask), top_n)
            se_df_subset = se_df[mask][:n] 
            
            print(f"Plotting {len(se_df_subset)} {target_ct} enriched SEs")
                  
        pdf_path = f"{outdir}/{data_source}_{_safe(target_ct)}_{len(se_df_subset)}_SEs_{direction}_specific{specific}.pdf"

        with PdfPages(pdf_path) as pdf:
            for idx, row in se_df_subset.iterrows():
                gene = row['gene_name']
                is_specific = row['is_specific'] 
                event = "_".join(idx.split("_")[1:])
            
                r_vals = psi_corr_df.loc[idx, ct_cols].values.astype(float)
                fdr_vals = psi_fdr_df.loc[idx, ct_cols].values.astype(float)
                lower = psi_ci_lower_df.loc[idx, ct_cols].values.astype(float)
                upper = psi_ci_upper_df.loc[idx, ct_cols].values.astype(float)

                # expression correlation
                try:
                    r_corr_vals = expr_corr_df.loc[gene, ct_cols].values.astype(float)
                except KeyError:
                    r_corr_vals = np.zeros(len(ct_cols))

                exon_coords = row['exon_coords']
                exon_len = row['exon_len']

                TITLE_SIZE = 15
                LABEL_SIZE = 14
                TICK_SIZE = 14

                fig, axes = plt.subplots(1, 2, figsize=(20, 8), dpi=200)
                fig.suptitle(f'{target_ct}', fontsize=20)

                colors = ['#4A90D9' if c == target_ct else '#AED6F1' for c in ct_cols]
                x = np.arange(len(ct_cols))

                # --- left panel: PSI correlations ---
                yerr = np.array([
                    np.maximum(r_vals - lower, 0),
                    np.maximum(upper - r_vals, 0)
                ])
                axes[0].bar(x, r_vals, color=colors, edgecolor='black', linewidth=0.8,
                            yerr=yerr, capsize=4, error_kw=dict(elinewidth=1, ecolor='black'))

                star_ys = []  # track y-position of every plotted "*"
                
                for i, fdr in enumerate(fdr_vals):
                    if fdr < fdr_thresh:
                        if r_vals[i] >= 0:
                            y_pos = upper[i] + 0.01
                            va = 'bottom'
                        else:
                            y_pos = lower[i] - 0.06
                            va = 'top'
                        
                        marker = "*"
                        if ct_cols[i] == target_ct and is_specific:
                            marker = "**"
                            
                        axes[0].text(i, y_pos, marker, ha='center', va=va, fontsize=18, color='black')
                        star_ys.append(y_pos)

                target_idx = ct_cols.index(target_ct)
                axes[0].set_title(
                    f'{gene} {event} (Spearman r = {np.round(r_vals[target_idx], 2)})\n'
                    f'{exon_len} nts ({exon_coords})',
                    fontsize=TITLE_SIZE
                )
                axes[0].set_ylabel('Spearman r', fontsize=LABEL_SIZE)
                axes[0].axhline(0, color='black', linewidth=0.6)
                axes[0].set_xticks(x)
                axes[0].set_xticklabels(ct_cols, rotation=45, ha='right', fontsize=TICK_SIZE)
                axes[0].tick_params(axis='y', labelsize=TICK_SIZE)

                # base extent from bars/CIs, extended to include star positions + margin
                y_min = min(-0.05, np.nanmin(lower))
                y_max = max(0.1, np.nanmax(upper) + 0.25)
                if star_ys:
                    star_margin = 0.15
                    y_min = min(y_min, min(star_ys) - star_margin)
                    y_max = max(y_max, max(star_ys) + star_margin)
                axes[0].set_ylim(y_min, y_max)

                # --- right panel: gene expression correlations ---
                axes[1].set_title(f'{gene} gene expression', fontsize=TITLE_SIZE)
                axes[1].bar(x, r_corr_vals, color=colors, edgecolor='black', linewidth=1)
                axes[1].axhline(0, color='black', linewidth=0.6)
                axes[1].set_xticks(x)
                axes[1].set_xticklabels(ct_cols, rotation=45, ha='right', fontsize=TICK_SIZE)
                axes[1].tick_params(axis='y', labelsize=TICK_SIZE)

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)