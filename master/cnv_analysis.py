from tp53_brca.main import *

# TODO: calculate all outliers without a cis effect and recheck samples with high numbers of outliers

# TODO: look at top # of outlier samples => check if CNV or high impact outliers overlap with outliers w/out cis effect ^
#  => look at top genes, top pathways, look for correlations with cancer pathways => is something 'out there'

# TODO: redo analysis per cohort


db_path = 'human_genome.db'

germline = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/cnv_germline_high_confidence.parquet')
splicing = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/fraser_aggregated_outliers_variants.parquet')
outrider = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/outrider_or_variants_predisppadjust_cnv.parquet')
protrider = pd.read_parquet('/s/project/cancer_pred/MASTER/final_res/protrider_pr_variants_predisppadjust_cnv.parquet')

sample_id = 'random_id'

germline_mapped = map_symbols_to_gene_ids(germline, db_path)

mrna_cnv = existing_cnv(outrider, germline_mapped, 'mRNA')
prot_cnv = existing_cnv(protrider, germline_mapped, 'Protein')


