from main import *
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 5000)

df = pd.read_csv('aberrant_expression_outliers.tsv', sep='\t')
db_path = 'human_genome.db'

results = find_trans_drivers(df)
results = get_chromosomal_relationship(results, db_path)

summary = results.copy()
summary['Has_CNV_Driver'] = summary['driver_gene'].notna()

ex = summary.groupby('sample_id')['Has_CNV_Driver'].any().value_counts(normalize=True)

print(f"Percentage of samples where trans-effects co-occur with at least one CNV: {ex[True]*100:.2f}%")

unique_drivers = pd.DataFrame(results['driver_gene'].unique(), columns=['geneID_short'])

driver_pathway_summary = get_top_pathways(unique_drivers, db_path, n=15)

print("=== TOP PATHWAYS OF IDENTIFIED DRIVERS ===")
print(driver_pathway_summary.to_string(index=False))

trans_chrom_hits = results[results['is_trans_chromosomal'] == True]

print("=== TOP CROSS-CHROMOSOME INTERACTIONS ===")
print(trans_chrom_hits[['sample_id', 'driver_gene', 'driver_chrom', 'target_gene', 'target_chrom']].head(10))
