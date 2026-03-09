import pandas as pd
import numpy as np
import pyarrow as par

'''
gtex_raw = pd.read_csv(
    '/home/emma/Downloads/GTEx_all_tissues_exported_counts.tsv',
    sep='\t',
    nrows=100000
)


gtex_raw.to_parquet('gtex.parquet')


root_cause = pd.read_csv('/home/emma/Downloads/results_cyclic_root_cause_analyis.tsv', sep='\t')

root_cause.to_parquet('root_cause.parquet')


root_cause = pd.read_parquet('root_cause.parquet')

numeric_df = root_cause.select_dtypes(include=[np.number])

matrix = numeric_df.fillna(0).values

rows, cols = np.where(matrix > 0.5)


hits = pd.DataFrame({
    'sampleID': root_cause['sampleID'].values[rows],
    'geneID': numeric_df.columns[cols],
    'probability': matrix[rows, cols]
})

hits = hits.sort_values(by='probability', ascending=False)


hits['clean_id'] = hits['geneID'].str.split('.').str[0]

hits.to_parquet('hits.parquet')


gtex = pd.read_parquet('gtex.parquet')
hits = pd.read_parquet('hits.parquet')

db_path = '/home/emma/PycharmProjects/thesis/gene_atlas/human_genome.db'

conn = sqlite3.connect(db_path)


functional_genes_query = """
    SELECT DISTINCT ensembl_gene_id 
    FROM genes 
    WHERE ensembl_gene_id IN (SELECT ensembl_id FROM pathways)
    OR ensembl_gene_id IN (SELECT gene1 FROM interactions)
"""

functional_genes = pd.read_sql(functional_genes_query, conn)['ensembl_gene_id'].tolist()
conn.close()


real_hits = hits[hits['clean_id'].isin(functional_genes)].copy()

merged_hits = pd.merge(
    real_hits,
    gtex,
    on=['sampleID', 'geneID'],
    how='inner'
)

merged_hits.to_parquet('merged.parquet')
'''
