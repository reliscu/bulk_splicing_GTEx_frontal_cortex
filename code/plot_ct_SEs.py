import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch

def _safe(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name))

def plot_ctype_SEs(
    ctype_specific_SEs,
    psi_corr_df,
    psi_fdr_df,
    expr_corr_df,
    SE_coords_df,
    psi_rank_centered,
    ct_ranks_dict,
    outdir,
    psi_ci_lower_df,   
    psi_ci_upper_df,
    fdr_thresh=0.05,
    ascending=False,
    n_boot=1000
):
    ct_cols = [c for c in psi_corr_df.columns if c != 'Gene']
    os.makedirs(outdir, exist_ok=True)
    psi_rank_np = np.array(psi_rank_centered, dtype=float)  # convert once

    for target_ct, SE_df in ctype_specific_SEs.items():
        if len(SE_df) == 0:
            print(f"Skipping {target_ct} — no specific SEs")
            continue

        pdf_path = f"{outdir}/{_safe(target_ct)}_specific_SEs_ascending{ascending}.pdf"
        print(f"Plotting {len(SE_df)} SEs for {target_ct} -> {pdf_path}")

        with PdfPages(pdf_path) as pdf:
            for idx, row in SE_df.iterrows():
                gene = row['Gene']
                SE = "_".join(idx.split("_")[1:])

                SE_pos = int(psi_corr_df.index.get_loc(idx))
                r_vals   = psi_corr_df.loc[idx, ct_cols].values.astype(float)
                fdr_vals = psi_fdr_df.loc[idx, ct_cols].values.astype(float)
                lower = psi_ci_lower_df.loc[idx, ct_cols].values.astype(float)
                upper = psi_ci_upper_df.loc[idx, ct_cols].values.astype(float)

                # expression correlation
                try:
                    r_corr_vals = expr_corr_df.loc[gene].values.astype(float)
                except KeyError:
                    r_corr_vals = np.zeros(len(ct_cols))

                # SE coordinates
                start_coord = int(SE_coords_df.loc[idx]['SE_start'])
                end_coord   = int(SE_coords_df.loc[idx]['SE_end'])
                SE_len    = end_coord - start_coord + 1

                fig, axes = plt.subplots(1, 2, figsize=(17, 6.5), dpi=200)
                fig.suptitle(f'{target_ct}', fontsize=16)

                colors = ['salmon' if c == target_ct else 'steelblue' for c in ct_cols]
                x = np.arange(len(ct_cols))

                # --- left panel: PSI correlations ---
                yerr = np.array([
                    np.maximum(r_vals - lower, 0),
                    np.maximum(upper - r_vals, 0)
                ])
                axes[0].bar(x, r_vals, color=colors, yerr=yerr, capsize=4,
                            error_kw=dict(elinewidth=1, ecolor='black'))

                for i, fdr in enumerate(fdr_vals):
                    if fdr < fdr_thresh:
                        if r_vals[i] >= 0:
                            y_pos = upper[i] + 0.02
                            va = 'bottom'
                        else:
                            y_pos = lower[i] - 0.02
                            va = 'top'
                        axes[0].text(i, y_pos, '*', ha='center', va=va,
                                     fontsize=12, color='black')

                target_idx = ct_cols.index(target_ct)
                axes[0].set_title(
                    f'{gene} {SE} SE (Spearman r = {np.round(r_vals[target_idx], 2)})\n'
                    f'{SE_len} nts ({start_coord}-{end_coord})'
                )
                axes[0].set_ylabel('Spearman r')
                axes[0].axhline(0, color='black', linewidth=0.6)
                axes[0].set_xticks(x)
                axes[0].set_xticklabels(ct_cols, rotation=45, ha='right', fontsize=9)
                axes[0].set_ylim(
                    min(-0.1, np.nanmin(lower) - 0.05),
                    max(0.1, np.nanmax(upper) + 0.1)
                )

                # --- right panel: gene expression correlations ---
                axes[1].set_title(f'{gene} gene expression')
                axes[1].bar(x, r_corr_vals, color=colors)
                axes[1].set_ylabel('Spearman r')
                axes[1].axhline(0, color='black', linewidth=0.6)
                axes[1].set_xticks(x)
                axes[1].set_xticklabels(ct_cols, rotation=45, ha='right', fontsize=9)

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)