import pandas as pd
import pyarrow as par

gtex_raw = pd.read_csv('/home/emma/Downloads/GTEx_all_tissues_exported_counts.tsv', sep='\t')

gtex_raw.to_parquet('gtex.parquet')

#print(gtex.head())