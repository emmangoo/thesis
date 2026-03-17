import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


root_cause = pd.read_parquet('root_cause.parquet')

print('root cause df: \n', root_cause.head())


numeric_df = root_cause.select_dtypes(include=[np.number])


matrix = numeric_df.fillna(0).values

all_values = matrix.flatten()
global_stats = pd.Series(all_values).describe(percentiles=[i/100 for i in range(90,101)])

print("Global Statistics across all cells:")
print(global_stats)

percentile_curve = global_stats.drop(['count', 'mean', 'std', 'min', 'max', '50%'])
percentile_curve.plot(kind='line', figsize=(10, 5))
plt.title("Top 10% Percentiles")
plt.grid(True, alpha=0.3)
plt.show()


top = pd.Series(all_values).describe(percentiles=[0.9, 0.95, 0.99, 0.999, 0.9999])
print(top)


'''
rows, cols = np.where(matrix > 0.5) # TODO: what number actually makes sense here


hits = pd.DataFrame({
    'sampleID': root_cause['sampleID'].values[rows],
    'geneID': numeric_df.columns[cols],
    'probability': matrix[rows, cols]
})

hits = hits.sort_values(by='probability', ascending=False)


hits['clean_id'] = hits['geneID'].str.split('.').str[0]

print(hits.head(), '\n hits shape: ', hits.shape)
'''