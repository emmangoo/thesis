from tp53_brca.main import *

# TODO: calculate all outliers without a cis effect and recheck samples with high numbers of outliers

# TODO: look at top # of outlier samples => check if CNV or high impact outliers overlap with outliers w/out cis effect ^
#  => look at top genes, top pathways, look for correlations with cancer pathways => is something 'out there'

# TODO: redo analysis per cohort


db_path = 'human_genome.db'

germline = pd.read_parquet('')
splicing = pd.read_parquet('')
outrider = pd.read_parquet('')
protrider = pd.read_parquet('')

sample_id = 'random_id'

germline_mapped = map_symbols_to_gene_ids(germline, db_path)

mrna_cnv = existing_cnv(outrider, germline_mapped, 'mRNA')
prot_cnv = existing_cnv(protrider, germline_mapped, 'Protein')

mrna_trans = mrna_cnv[mrna_cnv['Mechanism'] == 'trans effect'].copy()
prot_trans = prot_cnv[prot_cnv['Mechanism'] == 'trans effect'].copy()

def get_trans_burden(df, threshold=80):
    counts = df.groupby(sample_id).size().reset_index(name='trans_count')
    counts['is_hyper_trans'] = counts['trans_count'] >= threshold
    return counts

mrna_trans_burden = get_trans_burden(mrna_trans)
prot_trans_burden = get_trans_burden(prot_trans)

burden_ids_mrna = mrna_trans_burden[mrna_trans_burden['is_trans_outlier'] == True][sample_id].tolist()
burden_ids_prot = prot_trans_burden[prot_trans_burden['is_trans_outlier'] == True][sample_id].tolist()

print(f"number of samples identified as 'super outliers' (mrna): {len(burden_ids_mrna)}")
print(f"number of samples identified as 'super outliers' (prot): {len(burden_ids_prot)}")


def plot_entire_cohort_distribution(df):
    unique_samples = df.drop_duplicates(sample_id)
    total_patients = len(unique_samples)
    dist = unique_samples['Oncotree Code'].value_counts()
    dist_pct = unique_samples['Oncotree Code'].value_counts(normalize=True) * 100

    plt.figure(figsize=(12, 40))
    sns.countplot(data=unique_samples, y='Oncotree Code',
                  order=dist.index, palette='viridis')

    plt.title(f"Full Cohort Diagnosis Distribution (N = {total_patients})",
              fontsize=16, fontweight='bold')
    plt.xlabel("Number of Samples")
    plt.ylabel("Oncotree Code")
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.yticks(fontsize=8)

    plt.tight_layout()
    plt.savefig('cohort_distribution.png')
    plt.show()

    summary = pd.DataFrame({
        'Count': dist,
        'Percentage (%)': dist_pct
    })
    print("\n--- COHORT SUMMARY TABLE ---")
    print(summary)

    return unique_samples

cohort_df = plot_entire_cohort_distribution(mrna_cnv)