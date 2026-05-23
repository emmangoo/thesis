'''
Here, I am updating the database with BioMart chromosome/ scaffold names and karyotype bands,
as well as COSMIC data:
- https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/breakpoints
- https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/completecna
- https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/mutantcensus
- https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/completegeneexpression
- https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/hallmarksofcancer
- https://cancer.sanger.ac.uk/cosmic/download/cosmic/v103/structuralvariants
'''
import pandas as pd
import sqlite3


db_path = '/home/emma/PycharmProjects/thesis/gene_atlas/human_genome.db'
conn = sqlite3.connect(db_path)


def update_database():
    print("Updating BioMart data...")
    mart_path = '/home/emma/Downloads/database_extension/mart_export.txt'
    mart = pd.read_csv(mart_path, sep='\t', dtype={'Chromosome/scaffold name': str})

    mart_subset = mart[['Gene stable ID', 'Chromosome/scaffold name', 'Karyotype band']].drop_duplicates(
        'Gene stable ID')

    genes_df = pd.read_sql("SELECT * FROM genes", conn)
    updated_genes = pd.merge(genes_df, mart_subset,
                             left_on='ensembl_gene_id',
                             right_on='Gene stable ID',
                             how='left')

    updated_genes.to_sql('genes', conn, if_exists='replace', index=False)
    print("Genes table updated with BioMart info.")

    cosmic_files = {
        'cosmic_cna': '/home/emma/Downloads/database_extension/Cosmic_CompleteCNA_Tsv_v103_GRCh38/Cosmic_CompleteCNA_v103_GRCh38.tsv',
        'cosmic_mutants': '/home/emma/Downloads/database_extension/Cosmic_MutantCensus_Tsv_v103_GRCh38/Cosmic_MutantCensus_v103_GRCh38.tsv',
        'cosmic_hallmarks': '/home/emma/Downloads/database_extension/Cosmic_CancerGeneCensusHallmarksOfCancer_Tsv_v103_GRCh38/Cosmic_CancerGeneCensusHallmarksOfCancer_v103_GRCh38.tsv',
        'cosmic_struct_variants': '/home/emma/Downloads/database_extension/Cosmic_StructuralVariants_Tsv_v103_GRCh38/Cosmic_StructuralVariants_v103_GRCh38.tsv',
        'cosmic_breakpoints': '/home/emma/Downloads/database_extension/Cosmic_Breakpoints_Tsv_v103_GRCh38/Cosmic_Breakpoints_v103_GRCh38.tsv'
    }

    for table_name, path in cosmic_files.items():
        print(f"Processing {table_name}...")
        df = pd.read_csv(path, sep='\t', low_memory=False)
        df.columns = [c.lower().replace(' ', '_') for c in df.columns]
        df.to_sql(table_name, conn, if_exists='replace', index=False)

        cursor = conn.cursor()
        cols = df.columns.tolist()
        idx_col = 'gene_symbol' if 'gene_symbol' in cols else ('gene_name' if 'gene_name' in cols else None)

        if idx_col:
            cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_gene ON {table_name}({idx_col})")
            print(f"  Indexed {table_name} on {idx_col}")

    conn.commit()


def stream_large_csv(file_path, table_name, chunk_size=100000):
    print(f"Streaming {file_path} into table '{table_name}'...")

    reader = pd.read_csv(file_path, sep='\t', chunksize=chunk_size, low_memory=False, engine='c')

    for i, chunk in enumerate(reader):
        chunk.columns = [c.lower().replace(' ', '_') for c in chunk.columns]
        mode = 'replace' if i == 0 else 'append'
        chunk.to_sql(table_name, conn, if_exists=mode, index=False)

        if i % 10 == 0:
            print(f"  Processed {i * chunk_size} rows...")

    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({table_name})")
    cols = [info[1] for info in cursor.fetchall()]

    idx_col = 'gene_symbol' if 'gene_symbol' in cols else ('gene_name' if 'gene_name' in cols else None)
    if idx_col:
        print(f"Creating index for {table_name} on {idx_col}...")
        cursor.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_gene ON {table_name}({idx_col})")

    conn.commit()
    print(f"Finished loading {table_name}.")


update_database()

gene_exp_path = '/home/emma/Downloads/database_extension/Cosmic_CompleteGeneExpression_Tsv_v103_GRCh38/Cosmic_CompleteGeneExpression_v103_GRCh38.tsv'
stream_large_csv(gene_exp_path, 'cosmic_gene_expression')

conn.close()
print("Database fully updated.")