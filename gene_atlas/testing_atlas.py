import sqlite3
import pandas as pd

pd.set_option('display.max_colwidth', None)

def get_gene_info(db_path, symb, top_int):
    conn = sqlite3.connect(db_path)

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


def get_gene_info_by_ensg(db_path, raw_ensg_id, top_int):

    clean_id = raw_ensg_id.split('.')[0]

    conn = sqlite3.connect(db_path)

    gene_query = f"SELECT * FROM genes WHERE ensembl_gene_id = '{clean_id}'"
    gene_df = pd.read_sql_query(gene_query, conn)

    if gene_df.empty:
        print(f"Ensembl ID {clean_id} not found in database.")
        conn.close()
        return

    symb = gene_df.iloc[0]['symbol']


    pathway_query = f"SELECT pathway_name, pathway_id FROM pathways WHERE ensembl_id = '{clean_id}'"
    pathways_df = pd.read_sql_query(pathway_query, conn)


    interactor_query = f"""
    SELECT g.symbol as partner_symbol, i.combined_score 
    FROM interactions i 
    JOIN genes g ON i.gene2 = g.ensembl_gene_id 
    WHERE i.gene1 = '{clean_id}' 
    ORDER BY i.combined_score DESC 
    LIMIT {top_int}
    """
    interactors_df = pd.read_sql_query(interactor_query, conn)

    conn.close()

    print(f"=== GENE INFORMATION: {symb} ({clean_id}) ===")
    print(f"Name: {gene_df.iloc[0]['name']}")
    print(f"Notes: {gene_df.iloc[0]['notes']}")
    print("\n--- Top Pathways ---")
    if not pathways_df.empty:
        print(pathways_df['pathway_name'].head(10).to_string(index=False))
    else:
        print("No pathways found.")

    print(f"\n--- Top {top_int} Interacting Partners ---")
    if not interactors_df.empty:
        print(interactors_df.to_string(index=False))
    else:
        print("No interactions found.")


#get_gene_info('AVP', 25)

#get_gene_info_by_ensg('ENSG00000235174.1_4', 25)