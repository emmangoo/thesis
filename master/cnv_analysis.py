from main import *

# done: calculate all outliers without a cis effect and recheck samples with high numbers of outliers

# done: look at top # of outlier samples => check if CNV or high impact outliers overlap with outliers w/out cis effect ^
#  => look at top genes, top pathways, look for correlations with cancer pathways => is something 'out there'

# TODO: redo analysis per cohort


db_path = 'human_genome.db'

germline = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/cnv_germline_high_confidence.parquet')
splicing = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/fraser_aggregated_outliers_variants.parquet')
outrider = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/outrider_or_variants_predisppadjust_cnv.parquet')
protrider = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/protrider_pr_variants_predisppadjust_cnv.parquet')

sample_id = 'random_id'

germline_mapped = map_symbols_to_gene_ids(germline, db_path)

mrna_cnv = existing_cnv(outrider, germline_mapped, 'mRNA')
prot_cnv = existing_cnv(protrider, germline_mapped, 'Protein')

mrna_trans = mrna_cnv[mrna_cnv['Mechanism'] == 'trans effect'].copy()
prot_trans = prot_cnv[prot_cnv['Mechanism'] == 'trans effect'].copy()

def get_top_trans_samples(trans_df, n=10):
    counts = trans_df.groupby(sample_id).size().sort_values(ascending=False)
    top_samples = counts.head(n).reset_index(name='trans_outlier_count')
    
    #print(f"Top {n} samples by number of outliers:")
    #print(top_samples)
    return top_samples

top_mrna_samples = get_top_trans_samples(mrna_trans)
top_prot_samples = get_top_trans_samples(prot_trans)


def profile_top_samples(trans_df, top_samples_df, db_path, name="mRNA"):
    target_ids = top_samples_df[sample_id].tolist()
    top_data = trans_df[trans_df[sample_id].isin(target_ids)]
    
    top_genes = (top_data.groupby('geneID_short')
                 .size()
                 .sort_values(ascending=False)
                 .head(20)
                 .reset_index(name='frequency_in_top_samples'))
    
    genes = pd.DataFrame(top_data['geneID_short'].unique(), columns=['geneID_short'])
    top_pathways = get_top_pathways(genes, db_path, n=15)
    
    print(f"\n--- {name} profile: top 10 samples ---")
    print(f"total unique genes: {len(genes)}")
    print("\ntop pathways:")
    print(top_pathways.to_string(index=False))
    
    return top_genes, top_pathways

mrna_genes, mrna_pathways = profile_top_samples(mrna_trans, top_mrna_samples, db_path, "mRNA")
prot_genes, prot_pathways = profile_top_samples(prot_trans, top_prot_samples, db_path, "Protein")


def find_trans_drivers_new(full_df):
    df = full_df[full_df['Mechanism'] == 'trans effect'][['random_id', 'geneID_short', 'Gene', 'zScore']]
    df.columns = [sample_id, 'target_id', 'target_symbol', 'target_zscore']

    driver = full_df[full_df['Mechanism'] == 'cis effect'][['random_id', 'geneID_short', 'Gene', 'combined_cnv', 'IMPACT_snv']]
    driver.columns = [sample_id, 'driver_id', 'driver_symbol', 'driver_CNV', 'driver_SNV']

    links = pd.merge(df, driver, on=sample_id, how='inner')
    links = links[links['target_id'] != links['driver_id']]

    return links


mrna_drivers = find_trans_drivers_new(mrna_cnv)
summary = (mrna_drivers.groupby(['driver_id', 'driver_symbol'])
               .agg(
                   total_targets_affected=('target_id', 'count'),
                   unique_samples=(sample_id, 'nunique')
               )
               .sort_values('unique_samples', ascending=False))
    
print("\n--- top 10 drivers ---")
print(summary.head(10))

summary_ready = summary.reset_index().rename(columns={'driver_id': 'geneID_short'})

pathw = get_top_pathways(summary_ready, db_path, n=15)
print('--- top 10 pathways of drivers ---')
print(pathw.to_string(index=False))
