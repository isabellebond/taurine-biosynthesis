import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

import matplotlib.pyplot as plt
import pandas as pd

def plot_one_barplot(
    df,
    col,
    label_col=None,
    title=None,
    ylabel=None,
    sort=False,
    figsize=(6, 8),
    pval_col=None,
    pval_threshold=0.05,
    color="steelblue"
):
    """
    Create a single horizontal barplot, with optional significance stars
    if a p-value column is provided.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe.
    col : str
        Column name for the barplot values.
    label_col : str, optional
        Column name for bar labels. If None, df.index is used.
    title : str, optional
        Title for the barplot.
    ylabel : str, optional
        Label for the y-axis.
    sort : bool, default False
        Whether to sort bars by values (ascending).
    figsize : tuple, default (6, 8)
        Size of the figure.
    pval_col : str, optional
        Column name containing p-values.
    pval_threshold : float, default 0.05
        Threshold for marking significance with a star.
    color : str, default "steelblue"
        Color of the bars.
    """

    # Labels
    labels = df[label_col] if label_col else df.index.astype(str)

    # Sort if requested
    if sort:
        df = df.sort_values(col, ascending=True)
        labels = df[label_col] if label_col else df.index.astype(str)

    # X limits
    x_min = df[col].min()
    x_max = df[col].max()
    padding = 0.05 * (x_max - x_min) if x_max != x_min else 0.1
    x_min -= padding
    x_max += padding

    fig, ax = plt.subplots(figsize=figsize)

    # Barplot
    bars = ax.barh(labels, df[col], color=color)
    ax.set_title(title if title else col)
    ax.set_ylabel(ylabel if ylabel else "")
    ax.set_xlim(x_min, x_max)

    # Clean look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(left=False, bottom=False)

    # Add stars for significance
    if pval_col is not None:
        for i, bar in enumerate(bars):
            if df.iloc[i][pval_col] < pval_threshold:
                width = bar.get_width()
                offset = 0.01 * (x_max - x_min)
                ax.annotate(
                    "*",
                    xy=(width + offset if width >= 0 else width - offset,
                        bar.get_y() + bar.get_height() / 2),
                    xytext=(0, 0),
                    textcoords="offset points",
                    ha="left" if width >= 0 else "right",
                    va="center",
                    fontsize=12,
                    color="black",
                    fontweight="bold"
                )

    plt.tight_layout()
    return fig, ax

def plot_two_barplots(
    df,
    col1,
    col2,
    label_col=None,
    left_title=None,
    right_title=None,
    ylabel=None,
    sort=False,
    figsize=(10, 8),
    pval_col1=None,
    pval_col2=None,
    pval_threshold=0.05
):
    """
    Create two horizontal barplots side by side sharing the same y-axis,
    with only the left y-axis labeled. Each bar is annotated with a star (*)
    if the corresponding p-value is below the threshold. The x-axis limits
    are matched based on the maximum value in col1 or col2.
    """

    # Labels
    labels = df[label_col] if label_col else df.index.astype(str)

    # Sort if requested
    if sort:
        df = df.sort_values(col1, ascending=True)  # ascending so largest at top in barh
        labels = df[label_col] if label_col else df.index.astype(str)

    # Determine max value for shared xlim
    x_max = max(df[col1].max(), df[col2].max()) * 1.05  # add 5% padding
    x_min = min(df[col1].min(), df[col2].min()) * 1.05  # add 5% padding
    if x_min > 0:
        x_min = 0
    elif x_max < 0:
        x_max = 0

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    # Left barplot
    bars1 = axes[0].barh(labels, df[col1], color="steelblue")
    axes[0].set_title(left_title if left_title else col1)
    axes[0].set_ylabel(ylabel if ylabel else "")

    axes[0].set_xlim(x_min, x_max)

    # Right barplot
    bars2 = axes[1].barh(labels, df[col2], color="indianred")
    axes[1].set_title(right_title if right_title else col2)
    axes[1].set_ylabel("")
    axes[1].tick_params(axis="y", labelleft=False)
    axes[1].set_xlim(x_min, x_max)

    # Hide spines and ticks for a clean look
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(left=False, bottom=False)

    # Add stars for significant values in left plot
    if pval_col1 is not None:
        for i, bar in enumerate(bars1):
            if df.iloc[i][pval_col1] < pval_threshold:
                width = bar.get_width()
                axes[0].annotate(
                    "*",
                    xy=(width, bar.get_y() + bar.get_height() / 2),
                    xytext=(3, 0),
                    textcoords="offset points",
                    ha="center",
                    va="center",
                    fontsize=12,
                    color="black",
                    fontweight="bold"
                )

    # Add stars for significant values in right plot
    if pval_col2 is not None:
        for i, bar in enumerate(bars2):
            if df.iloc[i][pval_col2] < pval_threshold:
                width = bar.get_width()
                axes[1].annotate(
                    "*",
                    xy=(width, bar.get_y() + bar.get_height() / 2),
                    xytext=(3, 0),
                    textcoords="offset points",
                    ha="center",
                    va="center",
                    fontsize=12,
                    color="black",
                    fontweight="bold"
                )

    plt.tight_layout()
    return fig, axes

for gene_set in ['GO_Molecular_Function_2025', 'GO_Cellular_Component_2025', 'GO_Biological_Process_2025', 'Reactome_Pathways_2024']:
    untreated = pd.read_csv('results/remove_431_TTN/gsea_WT_Untreated_vs_TTN_Untreated.filter_' + gene_set + '.txt', sep='\t')
    treated = pd.read_csv('results/remove_431_TTN/gsea_WT_Untreated_vs_TTN_Taurine.filter_' + gene_set + '.txt', sep='\t')

    all_terms = untreated.merge(treated, on='Term', suffixes=('_untreated', '_treated'))
    all_terms = all_terms.sort_values(by = 'FDR q-val_untreated', ascending=True)
    all_terms = all_terms.loc[all_terms['FDR q-val_untreated'] < 0.05]
    all_terms = all_terms.dropna(subset=['NES_untreated', 'NES_treated'])

    all_terms_untreated = all_terms.head(50)
    plot_one_barplot(
        all_terms_untreated,
        col="NES_untreated",
        label_col="Term",
        title="WT vs TTN Untreated",
        ylabel="NES",
        sort=True,
        figsize=(8, 12),
        pval_col="FDR q-val_untreated",
        color="steelblue"
    )
    plt.savefig(f'results/remove_431_TTN/barplot_WTvsTTN_untreated_{gene_set}.png', bbox_inches='tight')


    print(f"{len(all_terms)} significant terms in {gene_set},plotting only terms where treated is opposing direction or non-significant.")
    all_terms = all_terms.loc[(all_terms['NES_untreated'] * all_terms['NES_treated'] < 0) | (all_terms['FDR q-val_treated'] > 0.05)]
    print(f"{len(all_terms)} terms remain after filtering.")
    if len(all_terms) > 50:
        all_terms = all_terms.head(50)
    

    plot_two_barplots(
        all_terms,
        col1="NES_untreated",
        col2="NES_treated",
        pval_col1="FDR q-val_untreated",
        pval_col2="FDR q-val_treated",
        label_col="Term",
        left_title="WT vs TTN Untreated",
        right_title="WT vs TTN with Taurine",
        ylabel="NES",
        sort=True,
        figsize=(12, 8),
    )

    plt.savefig(f'results/remove_431_TTN/barplot_WTvsTTNtreated_and_untreated_{gene_set}.png', bbox_inches='tight')

    wt = pd.read_csv('results/remove_431_TTN/gsea_WT_Untreated_vs_WT_Taurine.filter_' + gene_set + '.txt', sep='\t')
    ttn = pd.read_csv('results/remove_431_TTN/gsea_TTN_Untreated_vs_TTN_Taurine.filter_' + gene_set + '.txt', sep='\t')

    all_terms = untreated.merge(treated, on='Term', suffixes=('_untreated', '_treated'))
    all_terms = all_terms.sort_values(by = 'FDR q-val_untreated', ascending=True)
    all_terms = all_terms.loc[all_terms['FDR q-val_untreated'] < 0.05]
    all_terms = all_terms.dropna(subset=['NES_untreated', 'NES_treated'])


    print(f"{len(all_terms)} significant terms in {gene_set},plotting only terms where treated is opposing direction or non-significant.")
    all_terms = all_terms.loc[(all_terms['FDR q-val_treated'] < 0.05) & (all_terms['NES_untreated'] * all_terms['NES_treated'] > 0)]
    print(f"{len(all_terms)} terms remain after filtering.")
    if len(all_terms) > 50:
        all_terms = all_terms.head(50)
    

    plot_two_barplots(
        all_terms,
        col1="NES_untreated",
        col2="NES_treated",
        pval_col1="FDR q-val_untreated",
        pval_col2="FDR q-val_treated",
        label_col="Term",
        left_title="WT vs WT with Taurine",
        right_title="TTN vs TTN with Taurine",
        ylabel="NES",
        sort=True,
        figsize=(12, 8),
    )

    plt.savefig(f'results/remove_431_TTN/barplot_taurinetreated_vs_untreated_{gene_set}.png', bbox_inches='tight')