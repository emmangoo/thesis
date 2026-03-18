import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sqlite3
import textwrap
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 5000)

"""
TO DO:
analyze underexpression and overexpression outliers separately:

zScore > 0 --> overexpression
zScore < 0 --> underexpression
You can also just focus on set of predisposition genes, which are all the genes that have a padjust_predisp_extended value.
sample_id is randomized sample_id
Oncotree Code is the cancer subtype. one can also check the pathways separately within each subtype.
Diag is a more general grouping of the samples, based on Tissue and Oncotree Code. This is another way to group samples, which might be useful too.
CNV shows the status of copy number variant
snv and indel columns have the info about rare small snv and indels.
"""

exp = pd.read_csv('aberrant_expression_outliers.tsv', sep='\t')

#print(exp.head(), '\n\n')

over = exp[exp['zScore'] > 0]
under = exp[exp['zScore'] < 0]

# just look at predisposition genes
over_predisp = over[over['padjust_predisp_extended'].notna()]
under_predisp = under[under['padjust_predisp_extended'].notna()]


db_path = '/home/emma/Desktop/thesis/human_genome.db'


def get_top_recurrent_genes(df, db_path, n=10):
    gene_counts = df.groupby('geneID_short').agg({
        'geneID_short': 'size',
        'zScore': 'mean'
    }).rename(columns={'sampleID': 'sample_count'})

    conn = sqlite3.connect(db_path)
    atlas_symbols = pd.read_sql("SELECT ensembl_gene_id as geneID_short, symbol, name FROM genes", conn)
    conn.close()

    result = gene_counts.merge(atlas_symbols, on='geneID_short')
    return result.sort_values('sample_count', ascending=False).head(n)


def get_top_pathways(df, db_path, n=10):
    all_ids = df['geneID_short'].unique().tolist()

    conn = sqlite3.connect(db_path)
    placeholders = ','.join(['?'] * len(all_ids))
    query = f"""
        SELECT pathway_name, COUNT(DISTINCT ensembl_id) as unique_outlier_genes
        FROM pathways 
        WHERE ensembl_id IN ({placeholders})
        GROUP BY pathway_name
        ORDER BY unique_outlier_genes DESC
        LIMIT {n}
    """
    pathway_summary = pd.read_sql(query, conn, params=all_ids)
    conn.close()

    return pathway_summary


def plot_pathway_comparison(data_dict, n=10):

    combined_list = []

    for label, df in data_dict.items():
        temp_df = df.head(n).copy()
        temp_df['Category'] = label
        combined_list.append(temp_df)

    plot_df = pd.concat(combined_list)

    plt.figure(figsize=(16, 8))
    sns.set_style("whitegrid")

    ax = sns.barplot(
        data=plot_df,
        y='pathway_name',
        x='unique_outlier_genes',
        hue='Category',
        palette='viridis'
    )

    plt.title(f'Top {n} Pathways Across Expression Categories', fontsize=15)
    plt.xlabel('Number of Unique Genes', fontsize=12)
    plt.ylabel('Pathway Name', fontsize=12)
    plt.legend(title='Expression Type', bbox_to_anchor=(1.02, 0), loc='lower right', borderaxespad=0.)
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.savefig('pathway_comparison.png')
    plt.show()


def plot_pathways_subtypes(full_df, db_path, top_n_pathways=5, top_n_subtypes=6):
    top_subtypes = (
        full_df['Oncotree Code']
        .value_counts()
        .head(top_n_subtypes)
        .index
        .tolist()
    )

    combined_list = []
    for ca_type in top_subtypes:
        subtype_raw = full_df[full_df['Oncotree Code'] == ca_type]
        pathway_summary = get_top_pathways(subtype_raw, db_path, n=top_n_pathways)
        pathway_summary['Subtype'] = ca_type
        combined_list.append(pathway_summary)

    plot_df = pd.concat(combined_list)

    orders = {
        s: (plot_df[plot_df['Subtype'] == s]
             .sort_values('unique_outlier_genes', ascending=False)['pathway_name']
             .tolist())
        for s in top_subtypes
    }

    g = sns.FacetGrid(
        data=plot_df,
        col="Subtype",
        col_wrap=3,
        sharey=False,
        height=6,
        aspect=2.5
    )

    def _barplot(data, **kwargs):
        subtype = data['Subtype'].iloc[0]
        order = orders[subtype]
        sns.barplot(
            data=data,
            x="unique_outlier_genes",
            y="pathway_name",
            order=order,
            palette="magma",
            **kwargs
        )

    g.map_dataframe(_barplot)

    for ax in g.axes.flat:
        ax.xaxis.set_tick_params(which='both', length=3)
        ax.tick_params(axis='x', labelbottom=True)
        title = ax.get_title().split("=")[1].strip()
        ax.set_title(f"Cancer: {title}", fontweight='bold')

    g.set_axis_labels("Unique Genes", "")
    plt.tight_layout()
    plt.savefig("pathway_subtypes.png")
    plt.show()



results_to_plot = {
    'All Expression Outliers': get_top_pathways(exp, db_path, n=10),
    'Overexpression (z score > 0)': get_top_pathways(over, db_path, n=10),
    'Underexpression (z score < 0)': get_top_pathways(under, db_path, n=10),
    'Overexpression only with predisposition genes': get_top_pathways(over_predisp, db_path, n=10),
    'Underexpression only with predisposition genes': get_top_pathways(under_predisp, db_path, n=10)
}

to_plot_small = {
    'Overexpression only \n with predisposition genes': get_top_pathways(over_predisp, db_path, n=10),
    'Underexpression only \n with predisposition genes': get_top_pathways(under_predisp, db_path, n=10)
}

plot_pathway_comparison(results_to_plot)
#plot_pathway_comparison(to_plot_small)

#plot_pathways_subtypes(exp, db_path, top_n_pathways=5, top_n_subtypes=6)



