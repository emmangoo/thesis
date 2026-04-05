from main import *


db_path = '/home/emma/Desktop/thesis/human_genome.db'

tsv = pd.read_csv('/home/emma/Downloads/aberrant_expression_outliers.tsv', sep='\t')

exp = tsv[tsv['padjust_predisp_extended'].notna()]


def analyse_pathway_cis_trans(df, db_path, keyword):
    conn = sqlite3.connect(db_path)
    query = f"SELECT DISTINCT ensembl_id FROM pathways WHERE pathway_name LIKE '%{keyword}%'"
    pathway_genes = pd.read_sql(query, conn)['ensembl_id'].tolist()
    conn.close()

    hits = df[df['geneID_short'].isin(pathway_genes)].copy()

    # idea: a cis effect will be a gene from the outlier set with a snv/ cnv value, because
    # that means the gene was locally "hit" => a trans effect will not have these values
    # because it was 'hit' from an upstream cascade

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

brca = analyse_pathway_cis_trans(exp, db_path, 'BRCA')
tp53 = analyse_pathway_cis_trans(exp, db_path, 'TP53')

'''
def check_molecular_chaos(story_df, title):
    summary = story_df.groupby(['Oncotree Code', 'Mechanism']).size().unstack(fill_value=0)

    # Calculate % Trans
    if 'trans effect' in summary.columns and 'cis effect' in summary.columns:
        summary['Pct_Trans'] = (summary['trans effect'] /
                                (summary['trans effect'] + summary['cis effect']) * 100)

    #print(f"--- {title} Analysis ---")
    #print(summary.sort_values('trans effect', ascending=False).head(10))
    return summary


brca_chaos = check_molecular_chaos(brca, "BRCA Pathway")
tp53_chaos = check_molecular_chaos(tp53, "TP53 Pathway")


def plot_chaos_stacked(summary_df, title, figsize):
    summary_df = summary_df.sort_values('Pct_Trans', ascending=False)

    ax = summary_df[['cis effect', 'trans effect']].plot(
        kind='bar',
        stacked=True,
        figsize=(figsize, 6),
        color=['#9fc8c8', '#298c8c']
    )

    plt.title(f'Cis and Trans Effects by Cancer Subtype: {title}', fontsize=15, fontweight='bold')
    plt.ylabel('Number of Outliers')
    plt.xlabel('Cancer Subtype (Oncotree Code)')
    plt.legend(['Cis Effect', 'Trans Effect'])
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


plot_chaos_stacked(brca_chaos, "BRCA Pathway", 18)
plot_chaos_stacked(tp53_chaos, "TP53 Pathway", 24)

def plot_zscore_distribution(story_df, title):
    plt.figure(figsize=(10, 6))
    sns.stripplot(
        data=story_df,
        x='Mechanism',
        y='zScore',
        jitter=True,
        alpha=0.6,
        palette={'cis effect': '#9fc8c8', 'trans effect': '#298c8c'}
    )
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    plt.title(f'Z Score by Cis or Trans Effect: {title}')
    plt.show()

plot_zscore_distribution(brca, "BRCA Pathway")
plot_zscore_distribution(tp53, "TP53 Pathway")
'''