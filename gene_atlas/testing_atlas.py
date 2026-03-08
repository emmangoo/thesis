import sqlite3
import pandas as pd

pd.set_option('display.max_colwidth', None)

def get_gene_info(symb, top_int):
    conn = sqlite3.connect('human_genome.db')

    # get general gene info
    gene_query = f"SELECT * FROM genes WHERE symbol = '{symb}'"
    gene_df = pd.read_sql_query(gene_query, conn)

    if gene_df.empty:
        print(f"Gene {symb} not found.")
        return

    ensg_id = gene_df.iloc[0]['ensembl_gene_id']

    # get pathways
    pathway_query = f"SELECT pathway_name, pathway_id FROM pathways WHERE ensembl_id = '{ensg_id}'"
    pathways_df = pd.read_sql_query(pathway_query, conn)

    # get top interactors
    interactor_query = f"SELECT g.symbol as partner_symbol, i.combined_score FROM interactions i JOIN genes g ON i.gene2 = g.ensembl_gene_id WHERE i.gene1 = '{ensg_id}' ORDER BY i.combined_score DESC LIMIT {top_int}"

    interactors_df = pd.read_sql_query(interactor_query, conn)

    conn.close()

    print(f"=== GENE INFORMATION: {symb} ===")
    print(f"Name: {gene_df.iloc[0]['name']}")
    print(f"Notes: {gene_df.iloc[0]['notes']}")
    print("\n--- Top Pathways ---")
    print(pathways_df['pathway_name'].head(10).to_string(index=False))
    print(f"\n--- Top {top_int} Interacting Partners ---")
    print(interactors_df.to_string(index=False))


get_gene_info('TP53', 25)