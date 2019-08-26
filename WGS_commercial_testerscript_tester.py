import unittest
from treatment_program_rcm_helper import *
from treatment_program_rcm import *

commercial_WGS_tester = commercial_WGS_tester('ugh','ugh','ugh',ignore=True)

break_down_mutation = lambda gene: commercial_WGS_tester.break_down_mutation(gene)
classify_COM_mutation =  lambda gene: commercial_WGS_tester.check_variant_commercial(break_down_mutation(gene))

class Test(unittest.TestCase):

	def test_regions_WGS(self):
		"""Tests if WGS regions are only those that we wanted to search and that all are represented"""
		WGS = pd.read_csv('WGS_aggregated_test_results',sep='\t')

		#Get regions detected by WGS test
		mutations = WGS[['mutation','drug']].values
		regions = [[mutation, break_down_mutation(mutation).gene_name, drug] for mutation, drug in mutations]

		#Now make sure no region that we found was not in list and build list of regions we found
		regions_found = {'INH':set(), 'FQ':set(),'SLIS':set(),'RIF':set()}
		for mutation, region, drug in regions:
			if(drug == 'ISONIAZID'):
				drug = 'INH'
			regions_look = drug_to_WGSregion[drug]
			self.assertTrue(region in regions_look, msg = 'If failed it failed with {} not in {} for {} and drug {} '.format(region, regions_look, mutation, drug))
			regions_found[drug].add(region)

		#Check we found all regions, we are not missing any regions
		for drug in ['INH','FQ','SLIS','RIF']:
			what_we_found = list(regions_found[drug])
			what_we_should_found = drug_to_WGSregion[drug]
			for region in what_we_should_found:
				self.assertTrue(region in what_we_found, msg='If failed we did not find {} region for drug {}'.format(region, drug))


	def test_regions_locations_commercial(self):
		"""Tests if commercial regions are only those that we wanted to search and that all are represented"""
		COM = pd.read_csv('commercial_aggregated_test_results',sep='\t')


		#Get regions detected by commercial test
		mutations = COM[['mutation','drug']].values
		regions = [[mutation, break_down_mutation(mutation).gene_name, drug, break_down_mutation(mutation).codon_location] for mutation, drug in mutations]

		#Now make sure no region that we found was not in list or not in position of interest and build list of regions we found
		regions_found = {'INH':set(), 'FQ':set(),'SLIS':set(),'RIF':set()}
		for mutation, region, drug, position in regions:
			if(drug == 'ISONIAZID'):
				drug = 'INH'
			regions_look = drug_to_COMregion[drug]
			self.assertTrue(region in regions_look, msg = 'If failed it failed with {} not in {} for {} and drug {} '.format(region, regions_look, mutation, drug))
			regions_found[drug].add(region)

			position_look = drug_to_COMposition[drug][region]
			self.assertTrue(int(position) in position_look, msg = 'If failed it failed with {} not in {} for {} drug and {} region'.format(position, position_look, drug, region))

		#Check we found all regions, we are not missing any regions
		for drug in ['INH','FQ','SLIS','RIF']:
			what_we_found = list(regions_found[drug])
			what_we_should_found = drug_to_COMregion[drug]
			for region in what_we_should_found:
				self.assertTrue(region in what_we_found, msg='If failed we did not find {} region for drug {}'.format(region, drug))

	def test_raw_duplicates_WGS(self):
		"""Test if have any duplicate entires in WGS raw results"""
		WGS = pd.read_csv('WGS_raw_test_results',sep='\t')
		duplication = WGS.duplicated().values
		self.assertTrue(True not in duplication, msg='If failed you have duplicate rows in WGS_raw_test_results')

	def test_raw_duplicates_COM(self):
		"""Test if have any duplicate entires in commercial raw results"""
		COM = pd.read_csv('commercial_raw_test_results',sep='\t')
		duplication = COM.duplicated().values
		self.assertTrue(True not in duplication, msg='If failed you have duplicate rows in COM_raw_test_results')

	# def test_phenotype_WGS(self):
	# 	"""Test to make sure all isolates have correct phenotype"""
	# 	WGS = pd.read_csv('WGS_raw_test_results',sep='\t')
	# 	strain_info = annotate(pd.read_csv('strain_info.tsv',sep='\t'))
	# 	for index, row in WGS.iterrows():
	# 		strain = row['strain']
	# 		drug = row['drug']
	# 		real_phenotype = strain_info[strain_info['strain'] == strain][drug].item()
	# 		if(float(row['resistant']) == 1.0):
	# 			self.assertTrue(real_phenotype == 1, msg = '{} was annotated resistant but not actually resistant annotated as {}'.format(strain, real_phenotype))
	# 		elif(float(row['susceptible']) == 1.0):
	# 			self.assertTrue(real_phenotype == 0, msg = '{} was annotated susceptible but not actually susceptible annotated as {}'.format(strain, real_phenotype))
	# 		else:
	# 			raise Exception("Isolate {} is neither resistant or susceptible to drug {} ".format(strain, drug))

	# def test_phenotype_COM(self):
	# 	"""Test to make sure all isolates have correct phenotype"""
	# 	COM = pd.read_csv('commercial_raw_test_results',sep='\t')
	# 	strain_info = annotate(pd.read_csv('strain_info.tsv',sep='\t'))
	# 	for index, row in COM.iterrows():
	# 		strain = row['strain']
	# 		drug = row['drug']
	# 		real_phenotype = strain_info[strain_info['strain'] == strain][drug].item()
	# 		if(float(row['resistant']) == 1.0):
	# 			self.assertTrue(real_phenotype == 1, msg = '{} was annotated resistant but not actually resistant annotated as {}'.format(strain, real_phenotype))
	# 		elif(float(row['susceptible']) == 1.0):
	# 			self.assertTrue(real_phenotype == 0, msg = '{} was annotated susceptible but not actually susceptible annotated as {}'.format(strain, real_phenotype))
	# 		else:
	# 			raise Exception("Isolate {} is neither resistant or susceptible to drug {} ".format(strain, drug))

	def test_classifier_COM_SLIS(self):
		"""Test commercial classifier to see if it can detect positive and negative for SLI commercial mutations"""

		for test,expected in [['SNP_N_1473247_C1402C_rrs',['synonymous','SLIS']],
							['SNP_N_1473247_C1402T_rrs',['asynonymous','SLIS']],
							['SNP_N_1473247_C1401C_rrs',['synonymous','SLIS']],
							['SNP_N_1473247_C1401T_rrs',['asynonymous','SLIS']],
							['SNP_N_1473247_C1403C_rrs',[False, False]],
							['SNP_N_1473247_C1400T_rrs',[False, False]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(10+2715332), ['asynonymous','SLIS']],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(10+2715332), ['synonymous','SLIS']],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(11+2715332), ['asynonymous','SLIS']],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(11+2715332), ['synonymous','SLIS']],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(12+2715332), ['asynonymous','SLIS']],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(12+2715332), ['synonymous','SLIS']],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(13+2715332), ['asynonymous','SLIS']],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(13+2715332), ['synonymous','SLIS']],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(14+2715332), ['asynonymous','SLIS']],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(14+2715332), ['synonymous','SLIS']],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(37+2715332), ['asynonymous','SLIS']],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(37+2715332), ['synonymous','SLIS']],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(38+2715332), [False, False]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(9+2715332), [False, False]],
							['LSP_CN_410008_GG1-1402CT_VGQA1401-1402VGA_rrs', ['asynonymous', 'SLIS']],
							['LSP_CN_410008_TTGG1401-1405G_VGQA1401-1402GA_rrs', ['asynonymous', 'SLIS']],
							['LSP_CN_410008_TTGG1401-1402G_VGQA1390-1405VGA_rrs', ['asynonymous', 'SLIS']],
							['LSP_CN_410008_TTGG1390-1400G_VGQA1390-1400VGA_rrs', [False, False]],
							['DEL_P_4243204_d1401CT_rrs', ['asynonymous', 'SLIS']],
							['DEL_P_4243204_d1399CT_rrs', [False, False]],
							['INS_CF_2154110_i2002T_668L_rrs',[False, False]],
							['INS_CF_2154110_i1401T_668L_rrs',['asynonymous', 'SLIS']]]:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify as {} reported as {}'.format(test,expected, result))

	def test_classifier_COM_FQ(self):
		"""Test commercial classifier to see if it can detect positive and negative for FLQ commercial mutations"""
		
		tests = []

		for i in range(500,542):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrB'.format(i),['asynonymous','FQ']])

		for i in range(500,542):
			tests.append(['SNP_CN_6737_A1615G_T{}T_gyrB'.format(i),['synonymous','FQ']])

		for i in range(0,490):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrB'.format(i),[False, False]])

		for i in range(543,1000):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrB'.format(i),[False, False]])

		for i in range(89,95):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrA'.format(i),['asynonymous','FQ']])

		for i in range(89,95):
			tests.append(['SNP_CN_6737_A1615G_T{}T_gyrA'.format(i),['synonymous','FQ']])

		for i in range(0,88):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrA'.format(i),[False, False]])

		for i in range(96,150):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrA'.format(i),[False, False]])

		for test,expected in tests:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify reported as {}'.format(test, result))

	def test_classifier_COM_RIF(self):
		"""Test commercial classifier to see if it can detect positive and negative for RIF commercial mutations"""
		
		tests = []

		for i in range(424,454):
			tests.append(['SNP_CN_761141_C1335G_H{}Q_rpoB'.format(i),['asynonymous','RIF']])

		for i in range(424,454):
			tests.append(['SNP_CN_761141_C1335G_H{}H_rpoB'.format(i),['synonymous','RIF']])

		for i in range(0,424):
			tests.append(['SNP_CN_6737_A1615G_T{}A_rpoB'.format(i),[False, False]])

		for i in range(455,1000):
			tests.append(['SNP_CN_6737_A1615G_T{}A_rpoB'.format(i),[False, False]])

		for test,expected in tests:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify reported as {}'.format(test, result))

	def test_classifier_COM_INH(self):
		"""Test commercial classifier to see if it can detect positive and negative for INH commercial mutations"""
		
		tests = []
		tests.append(['SNP_CN_761141_C1335G_S315T_katG',['asynonymous','ISONIAZID']])
		tests.append(['SNP_CN_761141_C1335G_S315S_katG',['synonymous','ISONIAZID']])
		tests.append(['SNP_CN_761141_C1335G_S314S_katG',[False, False]])
		tests.append(['SNP_CN_761141_C1335G_S316S_katG',[False, False]])
		tests.append(['SNP_CN_761141_C1335G_S400S_katG',[False, False]])
		tests.append(['SNP_P_1673432_T8C_promoter-fabG1-inhA',['asynonymous','ISONIAZID']])
		tests.append(['SNP_P_1673432_T8T_promoter-fabG1-inhA',['synonymous','ISONIAZID']])
		tests.append(['SNP_P_1673432_T15C_promoter-fabG1-inhA',['asynonymous','ISONIAZID']])
		tests.append(['SNP_P_1673432_T15T_promoter-fabG1-inhA',['synonymous','ISONIAZID']])
		tests.append(['SNP_P_1673432_T16C_promoter-fabG1-inhA',['asynonymous','ISONIAZID']])
		tests.append(['SNP_P_1673432_T16T_promoter-fabG1-inhA',['synonymous','ISONIAZID']])
		tests.append(['SNP_P_1673432_T17C_promoter-fabG1-inhA',[False, False]])
		tests.append(['SNP_P_1673432_T13T_promoter-fabG1-inhA',[False, False]])

		for test,expected in tests:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify as {} reported as {}'.format(test, expected, result))


if __name__ == '__main__':
	unittest.main()