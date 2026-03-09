import pandas as pd
import pyarrow as par

"""
gtex_raw = pd.read_csv(
    '/home/emma/Downloads/GTEx_all_tissues_exported_counts.tsv',
    sep='\t',
    nrows=100000
)


gtex_raw.to_parquet('gtex.parquet')
"""

root_cause = pd.read_csv('/home/emma/Downloads/results_cyclic_root_cause_analyis.tsv', sep='\t')

root_cause.to_parquet('root_cause.parquet')
