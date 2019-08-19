import pandas as pd

WGS_results = pd.read_csv('WGS_results',sep='\t')
commercial_results = pd.read_csv('commercial_results',sep='\t')

WGS_results['count'] = 1
WGS_results.groupby(['drug','mutation']).sum().sort_values('drug', ascending=False).to_csv('commercial',sep='\t')

commercial_results['count'] = 1
commercial_results.groupby(['drug','mutation']).sum().sort_values('drug', ascending=False).to_csv('WGS',sep='\t')


