import pandas as pd
import numpy as np
import sqlite3
from gene_atlas.testing_atlas import get_gene_info_by_ensg
import pyarrow as par
#pd.set_option('display.max_columns', None)
#pd.set_option('display.width', 5000)

gtex = pd.read_parquet('gtex.parquet')
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

top_100_hits = hits.nlargest(100, 'probability')

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


hits['clean_id'] = hits['geneID'].str.split('.').str[0]

real_hits = hits[hits['clean_id'].isin(functional_genes)].copy()

print(f"Found {len(real_hits)} hits with functional data (skipping pseudogenes/noise).")

for _, row in real_hits.head(10).iterrows():
    print("\n" + "="*50)
    get_gene_info_by_ensg(db_path, row['geneID'], top_int=5)
    print(f"Score: {row['probability']}")