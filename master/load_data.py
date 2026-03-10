import pandas as pd
import numpy as np
import sqlite3
import pyarrow as par
import pyarrow.parquet as pq
import duckdb

'''
# reading in just the first part:

gtex_raw = pd.read_csv(
    '/home/emma/Downloads/GTEx_all_tissues_exported_counts.tsv',
    sep='\t',
    nrows=100000
)


gtex_raw.to_parquet('gtex.parquet')
'''

'''
# reading in the full file (in chunks):

tsv_path = '/home/emma/Downloads/GTEx_all_tissues_exported_counts.tsv'
parquet_path = 'gtex_full.parquet'

chunk_size = 100000

reader = pd.read_csv(tsv_path, sep='\t', chunksize=chunk_size, low_memory=False)

first_chunk = next(reader)
table = par.Table.from_pandas(first_chunk)
writer = pq.ParquetWriter(parquet_path, table.schema)
writer.write_table(table)

print("started")

try:
    for i, chunk in enumerate(reader):
        table = par.Table.from_pandas(chunk)
        writer.write_table(table)

        if i % 10 == 0:
            print(f"Written {(i + 2) * chunk_size} rows...")
finally:
    writer.close()

print("done")



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


# this will only work with the smaller gtex file:
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

# read in the full parquet file

db_path = '/home/emma/PycharmProjects/thesis/gene_atlas/human_genome.db'
hits_path = 'hits.parquet'
gtex_path = 'gtex_full.parquet'
output_path = 'merged_functional.parquet'


con = duckdb.connect()
con.execute(f"INSTALL sqlite; LOAD sqlite;")
con.execute(f"ATTACH '{db_path}' AS atlas (TYPE SQLITE);")

print("Merging Parquet files and filtering by SQLite Atlas...")

query = f"""
    COPY (
        SELECT 
            h.sampleID, 
            h.geneID, 
            h.probability, 
            h.clean_id,
            g.read_count, 
            g.fpkm, 
            g.tissue
        FROM read_parquet('{hits_path}') h
        INNER JOIN read_parquet('{gtex_path}') g
            ON h.sampleID = g.sampleID 
            AND h.geneID = g.geneID
        WHERE h.clean_id IN (
            SELECT DISTINCT ensembl_gene_id FROM atlas.genes 
            WHERE ensembl_gene_id IN (SELECT ensembl_id FROM atlas.pathways)
            OR ensembl_gene_id IN (SELECT gene1 FROM atlas.interactions)
        )
    ) TO '{output_path}' (FORMAT 'PARQUET')
"""

con.execute(query)
con.close()
