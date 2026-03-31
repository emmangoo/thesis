from main import *
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 5000)
import seaborn as sns
import matplotlib.pyplot as plt

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


'''
output:

--- CHROMOSOMAL ANALYSIS ---
Percentage of Driver-Target pairs on different chromosomes: 94.73%
Percentage of samples where trans-effects co-occur with at least one CNV: 74.15%
=== TOP PATHWAYS OF IDENTIFIED DRIVERS ===
                                                             pathway_name  unique_outlier_genes
                                                 Neutrophil degranulation                   315
                                            Generic Transcription Pathway                   273
              Antigen processing: Ubiquitination & Proteasome degradation                   232
                                            mRNA Splicing - Major Pathway                   192
                                                              Neddylation                   184
            Major pathway of rRNA processing in the nucleolus and cytosol                   166
                                           Dengue Virus-Host Interactions                   166
                                          Separation of Sister Chromatids                   151
                                         Ub-specific processing proteases                   135
                              Regulation of expression of SLITs and ROBOs                   134
                                                        RAC1 GTPase cycle                   121
                                                     mRNA Polyadenylation                   120
                                             RHO GTPases Activate Formins                   108
Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)                   105
              SRP-dependent cotranslational protein targeting to membrane                   104
=== TOP CROSS-CHROMOSOME INTERACTIONS ===
    sample_id      driver_gene driver_chrom      target_gene target_chrom
0        2443  ENSG00000117475            1  ENSG00000089234           12
1        2443  ENSG00000116754            1  ENSG00000089234           12
2        2443  ENSG00000000457            1  ENSG00000089234           12
3        2443  ENSG00000066557            1  ENSG00000089234           12
4        2443  ENSG00000171865            2  ENSG00000089234           12
5        2443  ENSG00000170471           20  ENSG00000089234           12
7        2443  ENSG00000173653           11  ENSG00000089234           12
8        2443  ENSG00000171863            2  ENSG00000089234           12
9        2443  ENSG00000143443            1  ENSG00000089234           12
10       2443  ENSG00000204070           20  ENSG00000089234           12
'''