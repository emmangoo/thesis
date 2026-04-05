import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sqlite3
import re
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 5000)

exp = pd.read_csv('aberrant_expression_outliers.tsv', sep='\t')

#print(exp.head(), '\n\n')

over = exp[exp['zScore'] > 0]
under = exp[exp['zScore'] < 0]

# just look at predisposition genes
over_predisp = over[over['padjust_predisp_extended'].notna()]
under_predisp = under[under['padjust_predisp_extended'].notna()]


db_path = 'human_genome.db'


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


def get_significant_subtypes_pathways(full_df, db_path, top_n_pathways=10, min_records=40):
    counts = full_df['Oncotree Code'].value_counts()

    qualified_subtypes = counts[counts >= min_records].index.tolist()

    subtype_dict = {}

    for ca_type in qualified_subtypes:
        subtype_data = full_df[full_df['Oncotree Code'] == ca_type]

        pathway_summary = get_top_pathways(subtype_data, db_path, n=top_n_pathways)
        subtype_dict[ca_type] = pathway_summary

    return subtype_dict


def find_cancers_by_pathway_keyword(subtype_dict, keyword):
    matches = {}
    keyword = keyword.lower()

    for ca_type, df in subtype_dict.items():
        matched_rows = df[df['pathway_name'].str.contains(keyword, case=False, na=False)]

        if not matched_rows.empty:
            matches[ca_type] = matched_rows['pathway_name'].tolist()

    print(f"Cancers with pathways containing '{keyword}'")
    for ca, pathways in matches.items():
        print(f"{ca}: {len(pathways)} matches found.")

    return matches


def plot_pathway_heatmap(subtype_dict, top_n_pathways=30,
                         title=f"Top Pathways Across Cancer Subtypes", figsize_width=20):

    all_data = []
    for ca_type, df in subtype_dict.items():
        temp = df.copy()
        temp['Subtype'] = ca_type
        all_data.append(temp)

    long_df = pd.concat(all_data)


    pivot_df = long_df.pivot_table(
        index='pathway_name',
        columns='Subtype',
        values='unique_outlier_genes'
    ).fillna(0)


    top_pathways = pivot_df.sum(axis=1).sort_values(ascending=False).head(top_n_pathways).index
    filtered_pivot = pivot_df.loc[top_pathways]


    plt.figure(figsize=(figsize_width, 12))
    sns.heatmap(
        filtered_pivot,
        annot=True,  # Show the actual numbers in the cells
        fmt=".0f",
        cmap="YlGnBu",
        cbar_kws={'label': 'Number of Unique Outlier Genes'}
    )

    plt.title(title, fontsize=18, fontweight='bold')
    plt.ylabel("")
    plt.xlabel("Cancer Subtype (Oncotree Code)")
    plt.tight_layout()
    plt.savefig('pathway_heatmap.png')
    plt.show()

def cis_trans(hits):
    # overexpression zScore > 0, underexpression zScore < 0
    over = hits['zScore'] > 0
    under = hits['zScore'] < 0

    # underexpression: impact "high"
    snv_cis = (under) & (
            (hits['IMPACT_snv'] == 'HIGH') |
            (hits['IMPACT_indel'] == 'HIGH')
    )

    # overexpression: amplification or duplication
    # underexpression: deletion CNV (germline/ somatic/ heterozygous)
    cnv_over_cis = (over) & (hits['CNV'].str.contains('AMP|DUP', case=False, na=False))
    cnv_under_cis = (under) & (hits['CNV'].str.contains('DEL', case=False, na=False))

    cnv_cis = cnv_over_cis | cnv_under_cis

    hits['Mechanism'] = 'trans effect'
    hits.loc[snv_cis | cnv_cis, 'Mechanism'] = 'cis effect'

    return hits


def find_trans_drivers(df):
    trans = df[df['Mechanism'] == 'trans effect'][['sample_id', 'geneID_short', 'zScore']]
    trans.columns = ['sample_id', 'target_gene', 'target_zScore']

    cnv = df[df['CNV'] != 'No CNV'][['sample_id', 'geneID_short', 'CNV']]
    cnv.columns = ['sample_id', 'driver_gene', 'driver_CNV']

    # maybe use different method if too slow
    # left join because we want to keep the samples where no CNV outlier was found
    drivers = pd.merge(trans, cnv, on='sample_id', how='left')
    # for safety:
    drivers = drivers[drivers['target_gene'] != drivers['driver_gene']]

    return drivers


def get_chromosomal_relationship(results_df, db_path):
    conn = sqlite3.connect(db_path)
    coords = pd.read_sql("SELECT ensembl_gene_id as geneID_short, location FROM genes", conn)
    conn.close()

    def extract_chrom(loc_str):
        if pd.isna(loc_str) or loc_str == '':
            return None
        match = re.match(r'^([0-9,X,Y,M,T]+)', str(loc_str))
        return match.group(1) if match else None

    coords['chrom'] = coords['location'].apply(extract_chrom)

    results_df = results_df.merge(coords[['geneID_short', 'chrom']],
                                  left_on='driver_gene', right_on='geneID_short', how='left')
    results_df = results_df.rename(columns={'chrom': 'driver_chrom'}).drop(columns=['geneID_short'])

    results_df = results_df.merge(coords[['geneID_short', 'chrom']],
                                  left_on='target_gene', right_on='geneID_short', how='left')
    results_df = results_df.rename(columns={'chrom': 'target_chrom'}).drop(columns=['geneID_short'])

    results_df['is_trans_chromosomal'] = results_df['driver_chrom'] != results_df['target_chrom']

    valid_pairs = results_df.dropna(subset=['driver_chrom', 'target_chrom'])
    if not valid_pairs.empty:
        trans_pc = valid_pairs['is_trans_chromosomal'].mean() * 100
        print(f"\n--- CHROMOSOMAL ANALYSIS ---")
        print(f"Percentage of Driver-Target pairs on different chromosomes: {trans_pc:.2f}%")
    return results_df


# check if driver gene and target gene interact with each other
def check_direct_trans_interaction(driver_symb, target_symb, db_path):
    conn = sqlite3.connect(db_path)
    query = f"""
        SELECT i.combined_score 
        FROM interactions i
        JOIN genes g1 ON i.gene1 = g1.ensembl_gene_id
        JOIN genes g2 ON i.gene2 = g2.ensembl_gene_id
        WHERE g1.symbol = '{driver_symb}' AND g2.symbol = '{target_symb}'
    """
    score = pd.read_sql(query, conn)
    conn.close()

    if not score.empty:
        return score.iloc[0]['combined_score']
    return None
