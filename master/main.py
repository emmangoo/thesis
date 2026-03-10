import pandas as pd
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
from gene_atlas.testing_atlas import get_gene_info_by_ensg
import pyarrow as par
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 5000)

# print(real_hits['clean_id'].nunique()) => 15490
# print(hits['clean_id'].nunique()) => 27012

merged = pd.read_parquet('merged_functional.parquet')

print(merged.head(), '\n', merged.shape, merged['sampleID'].nunique())

outliers = merged.groupby('sampleID').size().sort_values(ascending=False)


# top 5%
cutoff = int(len(outliers)*0.05)
outlier_samples = outliers.head(cutoff).index.tolist()

print(f"Top 5%: Samples with more than {outliers.iloc[cutoff]} outliers")

# find common pathways
def get_sample_pathway_profile(sample_id):
    sample_hits = merged[merged['sampleID'] == sample_id]['clean_id'].tolist()

    conn = sqlite3.connect('/home/emma/Desktop/thesis/human_genome.db')
    placeholders = ','.join(['?'] * len(sample_hits))
    query = f"SELECT pathway_name, COUNT(*) as count FROM pathways WHERE ensembl_id IN ({placeholders}) GROUP BY pathway_name ORDER BY count DESC"

    common_pathways = pd.read_sql(query, conn, params=sample_hits)
    conn.close()
    return common_pathways

top_sample = outlier_samples[0]
nr_outliers = outliers.get(top_sample, 0)

print(f"\nPathway Profile for Sample {top_sample}:\nThis sample has {nr_outliers} outliers\n")
print(get_sample_pathway_profile(top_sample).head(10))



outliers.hist(bins=50)
plt.title('Distribution of outliers')
plt.xlabel('Number of outlier genes')
plt.ylabel('Number of samples')
plt.show()



