"""
Building a custom gene atlas

    Sources:
    - HGNC: https://www.genenames.org/download/#!/#tocAnchor-1-1
    - Reactome mapping data: https://reactome.org/download-data
    - STRING-db: https://string-db.org/cgi/download?sessionId=bWzRzw6OBxQk

"""

import pandas as pd
import sqlite3


hgnc = pd.read_csv('/home/emma/Desktop/thesis/hgnc_data.txt', sep='\t')
reactome = pd.read_csv('/home/emma/Desktop/thesis/reactome_data.txt', sep='\t',

names=['ensembl_id', 'pathway_id', 'url', 'pathway_name', 'evidence', 'species'])

reactome = reactome[reactome['species'] == 'Homo sapiens']
protein_links = pd.read_csv('/home/emma/Desktop/thesis/protein_links_data.txt', sep=' ')


# aliases for the protein ids
aliases = pd.read_csv('/home/emma/Desktop/thesis/aliases.txt', sep='\t')
gene_aliases = aliases[aliases['alias'].str.startswith('ENSG', na=False)].copy()

mapping_dict = dict(zip(gene_aliases['#string_protein_id'], gene_aliases['alias']))

# clean STRING-db file to only include highly confident protein-protein interactions
high_conf_links = protein_links[protein_links['combined_score'] >= 700].copy()


high_conf_links['gene1'] = high_conf_links['protein1'].map(mapping_dict)
high_conf_links['gene2'] = high_conf_links['protein2'].map(mapping_dict)

high_conf_links.dropna(subset=['gene1', 'gene2'], inplace=True)

conn = sqlite3.connect('human_genome.db')

hgnc.to_sql('genes', conn, if_exists='replace', index=False)
reactome.to_sql('pathways', conn, if_exists='replace', index=False)
high_conf_links.to_sql('interactions', conn, if_exists='replace', index=False)


cursor = conn.cursor()

try:
    cursor.execute("ALTER TABLE genes ADD COLUMN notes TEXT")
except sqlite3.OperationalError:
    pass


cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene1 ON interactions(gene1)")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene2 ON interactions(gene2)")
cursor.execute("CREATE INDEX IF NOT EXISTS idx_ensembl ON genes(ensembl_gene_id)")

conn.commit()

conn.close()
