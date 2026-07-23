from main import *
from scipy.stats import mannwhitneyu, fisher_exact, linregress


db_path = 'human_genome.db'

germline = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/cnv_germline_high_confidence.parquet')
splicing = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/fraser_aggregated_outliers_variants.parquet')
outrider = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/outrider_or_variants_predisppadjust_cnv.parquet')
protrider = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/protrider_pr_variants_predisppadjust_cnv.parquet')

sample_id = 'random_id'

# initial stats
def plot_initial_stats(df):
    plt.figure(figsize=(12,5))

    plt.subplot(1,2,1)
    sns.countplot(data=df, x='Type', palette='Set2', hue='Type', legend=False)
    plt.title('Frequency of CNV types')
    plt.xticks(rotation=45)

    plt.subplot(1,2,2)
    cnvs_per_sample = df.groupby(sample_id).size()
    sns.histplot(cnvs_per_sample, bins=30, kde=True, color='blue')
    plt.xlim(right=1500)
    plt.title('Distribution of CNV per sample')
    plt.xlabel('Number of CNVs')

    plt.tight_layout()
    plt.savefig('CNV_distribution.png')
    plt.show()

    return cnvs_per_sample.describe()

def get_top_cnv_genes(df, percentile=0.95):
    counts = (df.groupby('Gene')
              .agg({sample_id: 'nunique'})
              .rename(columns={sample_id: 'sample_counts'})
              .sort_values('sample_counts', ascending=False))

    top = counts['sample_counts'].quantile(percentile, interpolation='nearest')
    top_genes = counts[counts['sample_counts'] >= top]

    top_genes = top_genes.reset_index().rename(columns={'Gene' : 'geneID_short'})

    print(f'Found {len(top_genes)} genes (threshold {top} samples)')
    return top_genes


#plot_initial_stats(germline)
top_cnv_genes = get_top_cnv_genes(germline)
hotspots = get_top_pathways(top_cnv_genes, db_path, n=15)

def analyse_cnv_samples(cnv_df, splicing_df, mrna_df, protein_df, threshold=80):
    # outliers only (maybe adjust values if needed)
    splicing_out = splicing_df[splicing_df['padjustGene'] < 0.05]
    mrna_out = mrna_df[mrna_df['padjust'] < 0.05]
    protein_out = protein_df[protein_df['padjust'] < 0.05]

    all_outliers = (pd.concat([
        splicing_out[[sample_id]],
        mrna_out[[sample_id]],
        protein_out[[sample_id]]
    ])
    .groupby(sample_id).size()
    .reset_index(name='outlier_count'))

    all_cnv = cnv_df.groupby(sample_id).size().reset_index(name='cnv_count')

    comp = pd.merge(all_outliers, all_cnv, on=sample_id, how='left').fillna(0)

    comp['group'] = 'normal'
    comp.loc[comp['outlier_count'] >= threshold, 'group'] = 'outlier'

    outlier_cnvs = comp[comp['group'] == 'outlier']['cnv_count']
    normal_cnvs = comp[comp['group'] == 'normal']['cnv_count']

    stat, p = mannwhitneyu(outlier_cnvs, normal_cnvs, alternative='greater')

    plt.figure(figsize=(14,6))

    plt.subplot(1, 2, 1)
    sns.regplot(data=comp, x='outlier_count', y='cnv_count',
                scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})
    plt.ylim(top=1000) 
    plt.title(f"Outliers vs CNVs\n(Total Samples: {len(comp)})")
    plt.xlabel("Total Outlier Count")
    plt.ylabel("CNV Count (germline)")

    plt.subplot(1, 2, 2)
    sns.boxplot(data=comp, x='group', y='cnv_count', palette='Set2')
    plt.ylim(top=1000)
    plt.title(f"CNV Burden (P-Value: {p:.2e})")
    plt.ylabel("CNV Count")

    plt.tight_layout()
    plt.savefig('mannwhitneyu.png')
    plt.show()

    return comp, p


#comp, p = analyse_cnv_samples(germline, splicing, outrider, protrider)


germline_mapped = map_symbols_to_gene_ids(germline, db_path)

mrna_cnv = existing_cnv(outrider, germline_mapped, 'mRNA')
prot_cnv = existing_cnv(protrider, germline_mapped, 'Protein')



def get_cnv_hotspots(cnv, percentile=0.95):

    hotspot_stats = (cnv.groupby(['geneID_short', 'Gene'])
                      .agg(
                          sample_count=(sample_id, 'nunique'),
                          total_cnv_count=(sample_id, 'count')
                      )
                      .sort_values('sample_count', ascending=False)
                      ).reset_index()

    threshold = hotspot_stats['sample_count'].quantile(percentile)
    top = hotspot_stats[hotspot_stats['sample_count'] >= threshold].copy()

    print('Top 10 CNV hotspot genes \n')
    print(top[['Gene', 'sample_count', 'total_cnv_count', 'geneID_short']].head(10).to_string(index=False))
        
    print(f"{len(top)} hotspot genes in >= {threshold} samples")
    
    return top


hotspot_genes = get_cnv_hotspots(germline_mapped)
hotspot_pathways = get_top_pathways(hotspot_genes, db_path, n=15)

print('\n top pathways of cnv hotspots:\n')
print(hotspot_pathways.to_string(index=False))


def hotspot_significance(mrna_df, cnv_df, hotspot_genes, threshold=80):
    sample_counts = (mrna_df.groupby(sample_id)
                     .size()
                     .reset_index(name='outlier_count'))
    sample_counts['is_hyper_outlier'] = sample_counts['outlier_count'] >= threshold
    
    hotspot_ids = set(hotspot_genes['geneID_short'])
    hotspot_samples = set(cnv_df[cnv_df['geneID_short'].isin(hotspot_ids)][sample_id])
    
    sample_counts['has_hotspot_hit'] = sample_counts[sample_id].isin(hotspot_samples)
    
    # create 2x2 table with hotspot cnv / no hotspot cnv X hyper outlier / normal sample
    table = pd.crosstab(sample_counts['is_hyper_outlier'], sample_counts['has_hotspot_hit'])
    odds, p = fisher_exact(table)
    
    print("----- Results of Fisher's exact test ------")
    print(table)
    print(f"\n odds: {odds:.2f}")
    print(f"\n p-value: {p:.2e}")
    
    return table, odds, p

#table, odds, p = hotspot_significance(outrider, germline_mapped, hotspot_genes)


def plot_cn_zscore(outlier_df, cnv_df, layer_name="mRNA"):

    merged = pd.merge(
        outlier_df[['random_id', 'geneID_short', 'zScore', 'padjust']],
        cnv_df[['random_id', 'geneID_short', 'CN', 'Type']],
        on=['random_id', 'geneID_short'],
        how='inner'
    )

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=merged, x='CN', y='zScore', hue='Type', palette='Set2')
    
    plt.axhline(0, color='black', linestyle='--', alpha=0.5)
    
    plt.title(f"{layer_name} Copy Number vs. z Score", fontweight='bold')
    plt.xlabel("Copy Number")
    plt.ylabel("z Score")
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f'zscores_copynumber_{layer_name}')
    plt.show()

    correlation = merged['CN'].corr(merged['zScore'])
    print(f"Pearson Correlation between {layer_name} CN and Z-Score: {correlation:.3f}")


    clean_df = merged.dropna(subset=['CN', 'zScore'])
    
    slope, intercept, r_value, p_value, std_err = linregress(clean_df['CN'], clean_df['zScore'])
    
    print(f"Cis-Effect Slope: {slope:.4f}")
    print(f"R-squared: {r_value**2:.4f}")
    
    return merged

mrna_dose = plot_cn_zscore(outrider, germline_mapped, "mRNA")
prot_dose = plot_cn_zscore(protrider, germline_mapped, "Protein")

