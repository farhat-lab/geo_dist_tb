import pandas as pd

lineage_snps = [i for i in open('lineage_snp','r').readlines()]


WGS_results = pd.read_csv('WGS_results',sep='\t')
commercial_results = pd.read_csv('commercial_results',sep='\t')

WGS_results['count'] = 1
WGS_results = WGS_results[~WGS_results['mutation'].isin(lineage_snps)]
WGS_results.groupby(['drug','mutation']).sum().sort_values('drug', ascending=False).to_csv('WGS',sep='\t')

commercial_results['count'] = 1
commercial_results = commercial_results[~commercial_results['mutation'].isin(lineage_snps)]
commercial_results.groupby(['drug','mutation']).sum().sort_values('drug', ascending=False).to_csv('commercial',sep='\t')


