import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 5000)

"""
TO DO:
analyze underexpression and overexpression outliers separately:

zScore > 0 --> overexpression
zScore < 0 --> underexpression
You can also just focus on set of predisposition genes, which are all the genes that have a padjust_predisp_extended value.
sample_id is randomized sample_id
Oncotree Code is the cancer subtype. one can also check the pathways separately within each subtype.
Diag is a more general grouping of the samples, based on Tissue and Oncotree Code. This is another way to group samples, which might be useful too.
CNV shows the status of copy number variant
snv and indel columns have the info about rare small snv and indels.
"""

exp = pd.read_csv('aberrant_expression_outliers.tsv', sep='\t')

print(exp.head(), '\n\n')

over = exp[exp['zScore'] > 0]
under = exp[exp['zScore'] < 0]

# just look at predisposition genes
over_predisp = over[over['padjust_predisp_extended'].notna()]

