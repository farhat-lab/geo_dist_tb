import unittest
from treatment_program_rcm_helper import *
from treatment_program_rcm import *

commercial_WGS_tester_full = commercial_WGS_tester('/home/lf61/lf61/mic_assemblies/46-annotate-vcfs-yasha/flatann2','strain_info.tsv','results_modified_unknown', 'lineage_snp', 'snps','panel.final.Cryptic_no_frameshift.txt')
commercial_WGS_tester = commercial_WGS_tester('ugh','ugh','ugh','ugh','ugh','ugh',ignore=True)
break_down_mutation = lambda gene: commercial_WGS_tester.break_down_mutation(gene)
classify_COM_mutation =  lambda gene: commercial_WGS_tester.check_variant_commercial(break_down_mutation(gene))
classify_WGS_mutation = lambda gene: commercial_WGS_tester.check_variant_WGS(gene)


class Test(unittest.TestCase):

	def test_table_10_regions_WGS(self):
		"""Tests if WGS regions in supplementary table 10 are only those that we wanted to search and that all are represented"""
		WGS = pd.read_csv('WGS_aggregated_test_results',sep='\t')
		WGS_AJRCCM = WGS[WGS['extra_annotation'].str.contains("Table-10-snp")]

		#Get regions detected by WGS test
		mutations = WGS_AJRCCM[['mutation','drug']].values
		regions = [[mutation, break_down_mutation(mutation).gene_name, drug] for mutation, drug in mutations]

		#Now make sure no region that we found was not in list and build list of regions we found
		regions_found = {'INH':set(), 'FLQ':set(),'SLIS':set(),'RIF':set(), 'STR':set(),'EMB':set(),'PZA':set()}
		for mutation, region, drug in regions:
			if(drug == 'ISONIAZID'):
				drug = 'INH'
			if(drug == 'FQ'):
				drug = 'FLQ'
			regions_look = drug_to_WGSregion[drug]
			self.assertTrue(region in regions_look, msg = 'If failed it failed with {} not in {} for {} and drug {} '.format(region, regions_look, mutation, drug))
			regions_found[drug].add(region)

		#Check we found all regions, we are not missing any regions
		for drug in ['INH','FLQ','SLIS','RIF','EMB','PZA','STR']:
			what_we_found = list(regions_found[drug])
			what_we_should_found = drug_to_WGSregion[drug]
			for region in what_we_should_found:
				self.assertTrue(region in what_we_found, msg='If failed we did not find {} region for drug {}'.format(region, drug))

	def test_AJRCCM_labels(self):
		"""Test to make sure all AJRCCM labels are labeled correctly/none are not labeled that should be labeled """
		WGS = pd.read_csv('WGS_aggregated_test_results',sep='\t')
		mutations = [[drug,break_down_mutation(mutation), annotation] for drug,mutation, annotation in WGS[['drug','mutation','extra_annotation']].values]
		processed_snps = []
		snps = [i.split('\t') for i in open('snps','r').readlines()]
		for drug, snp in snps:
			if(drug == 'KAN' or drug == 'AMK' or drug == 'CAP'):
				drug = 'SLIS'
			elif(drug == 'LEVO' or drug == 'MOXI' or drug == 'OFLX'):
				drug = 'FLQ'
			elif(drug != 'INH' and drug != 'RIF' and drug != 'PZA' and drug != 'STR' and drug != 'EMB'):
				drug = None
			if(drug):
				processed_snps.append([break_down_mutation(snp.rstrip()), drug])

		for drug_compare, mutation, annotation in mutations:
			found = False
			for snp, drug in processed_snps:
				if(snp.compare_variant_name_location(mutation) and drug == drug_compare):
					found = True
			if(not found):
				if(mutation.gene_name in ['rpoB','pncA', 'katG']):
					change = mutation.AA_change
					if('-' in change):
						one, two = change.split('-')
						if(len(one) != len(two) and (len(one)%3!= 0 or len(two)%3!=0)):
							found = True
					else:
						if(len(change) == 1 or len(change) > 2):
							found = True

			if(found):
				self.assertTrue(found == True and ('AJRCCM' in annotation) , msg='Found was {} annotation was {} for {} which is not a correct label'.format(found, annotation, mutation))
			else:
				self.assertTrue(not found and ('AJRCCM' not in annotation), msg='Found was {} annotation was {} for {} which is not a correct label'.format(found, annotation, mutation))


	def test_CRYPTIC_labels(self):
		"""Test to make sure all CRYPTIC labels are labeled correctly/none are not labeled that should be labeled """
		WGS = pd.read_csv('WGS_aggregated_test_results',sep='\t')
		mutations = [[drug,break_down_mutation(mutation), annotation] for drug,mutation, annotation in WGS[['drug','mutation','extra_annotation']].values]
		processed_snps = []
		processed_snps = [commercial_WGS_tester.cryptic_snp_parser(i.rstrip()) for i in open('panel.final.Cryptic_no_frameshift.txt','r').readlines()]
		for drug_compare, mutation, annotation in mutations:
			found = False
			for drug, snp in processed_snps:
				if(snp.compare_variant_name_location(mutation) and drug_compare in drug):
					found = True
			if(not found):
				if(mutation.gene_name in ['rpoB','pncA', 'katG']):
					change = mutation.AA_change
					if('-' in change):
						one, two = change.split('-')
						if(len(one) != len(two) and (len(one)%3!= 0 or len(two)%3!=0)):
							found = True
					else:
						if(len(change) == 1 or len(change) > 2):
							found = True

			if(found):
				self.assertTrue(found == True and ('CRYPTIC' in annotation) , msg='Found was {} annotation was {} for {} which is not a correct label'.format(found, annotation, mutation))
			else:
				self.assertTrue(not found and ('CRYPTIC' not in annotation), msg='Found was {} annotation was {} for {} which is not a correct label'.format(found, annotation, mutation))



	def test_table_10_labels(self):
		"""Test to make sure all table 10 labels are labeled correctly/none are not labeled that should be """
		WGS = pd.read_csv('WGS_aggregated_test_results',sep='\t')
		# WGS_AJRCCM = WGS[WGS['extra_annotation']]
		mutations = [[break_down_mutation(mutation), annotation] for mutation, annotation in WGS[['mutation','extra_annotation']].values]

		for mutation, annotation in mutations:
			found = False
			if(classify_WGS_mutation(mutation)):
				found = True
			if(found):
				self.assertTrue(found == True and ('Table-10-snp' in annotation) or ('IGNORE' in annotation), msg='Found was {} annotation was {} for {} which is not a correct label'.format(found, annotation, mutation))
			else:
				self.assertTrue(not found and ('Table-10-snp' not in annotation) or ('IGNORE' in annotation), msg='Found was {} annotation was {} for {} which is not a correct label'.format(found, annotation, mutation))

	def test_regions_locations_commercial(self):
		"""Tests if commercial regions are only those that we wanted to search and that all are represented"""
		COM = pd.read_csv('commercial_aggregated_test_results',sep='\t')


		#Get regions detected by commercial test
		mutations = COM[['mutation','drug']].values
		regions = [[mutation, break_down_mutation(mutation).gene_name, drug, break_down_mutation(mutation).codon_location] for mutation, drug in mutations]

		#Now make sure no region that we found was not in list or not in position of interest and build list of regions we found
		regions_found = {'INH':set(), 'FLQ':set(),'SLIS':set(),'RIF':set()}
		for mutation, region, drug, position in regions:
			if(drug == 'ISONIAZID'):
				drug = 'INH'
			if(drug == 'FQ'):
				drug = 'FLQ'
			regions_look = drug_to_COMregion[drug]
			self.assertTrue(region in regions_look, msg = 'If failed it failed with {} not in {} for {} and drug {} '.format(region, regions_look, mutation, drug))
			regions_found[drug].add(region)

			position_look = drug_to_COMposition[drug][region]
			self.assertTrue(int(position) in position_look, msg = 'If failed it failed with {} not in {} for {} drug and {} region'.format(position, position_look, drug, region))

		#Check we found all regions, we are not missing any regions
		for drug in ['INH','FLQ','SLIS','RIF']:
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

		for test,expected in [['SNP_N_1473247_C1402C_rrs',['synonymous',['SLIS']]],
							['SNP_N_1473247_C1402T_rrs',['asynonymous',['SLIS']]],
							['SNP_N_1473247_C1401C_rrs',['synonymous',['SLIS']]],
							['SNP_N_1473247_C1401T_rrs',['asynonymous',['SLIS']]],
							['SNP_N_1473247_C1403C_rrs',[False, [False]]],
							['SNP_N_1473247_C1400T_rrs',[False, [False]]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(10+2715332), ['asynonymous',['SLIS']]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(10+2715332), ['synonymous',['SLIS']]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(11+2715332), ['asynonymous',['SLIS']]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(11+2715332), ['synonymous',['SLIS']]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(12+2715332), ['asynonymous',['SLIS']]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(12+2715332), ['synonymous',['SLIS']]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(13+2715332), ['asynonymous',['SLIS']]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(13+2715332), ['synonymous',['SLIS']]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(14+2715332), ['asynonymous',['SLIS']]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(14+2715332), ['synonymous',['SLIS']]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(37+2715332), ['asynonymous',['SLIS']]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(37+2715332), ['synonymous',['SLIS']]],
							['SNP_I_{}_C103A_inter-eis-Rv2417c'.format(38+2715332), [False, [False]]],
							['SNP_I_{}_C103C_inter-eis-Rv2417c'.format(9+2715332), [False, [False]]],
							['LSP_CN_410008_GG1-1402CT_VGQA1401-1402VGA_rrs', ['asynonymous', ['SLIS']]],
							['LSP_CN_410008_TTGG1401-1405G_VGQA1401-1402GA_rrs', ['asynonymous', ['SLIS']]],
							['LSP_CN_410008_TTGG1401-1402G_VGQA1390-1405VGA_rrs', ['asynonymous', ['SLIS']]],
							['LSP_CN_410008_TTGG1390-1400G_VGQA1390-1400VGA_rrs', [False, [False]]],
							['DEL_P_4243204_d1401CT_rrs', ['asynonymous', ['SLIS']]],
							['DEL_P_4243204_d1399CT_rrs', [False, [False]]],
							['INS_CF_2154110_i2002T_668L_rrs',[False, [False]]],
							['INS_CF_2154110_i1401T_668L_rrs',['asynonymous', ['SLIS']]]]:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify as {} reported as {}'.format(test,expected, result))

	def test_classifier_COM_FLQ(self):
		"""Test commercial classifier to see if it can detect positive and negative for FLQ commercial mutations"""
		
		tests = []

		for i in range(500,542):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrB'.format(i),['asynonymous',['FLQ']]])

		for i in range(500,542):
			tests.append(['SNP_CN_6737_A1615G_T{}T_gyrB'.format(i),['synonymous',['FLQ']]])

		for i in range(0,490):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrB'.format(i),[False, [False]]])

		for i in range(543,1000):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrB'.format(i),[False, [False]]])

		for i in range(89,95):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrA'.format(i),['asynonymous',['FLQ']]])

		for i in range(89,95):
			tests.append(['SNP_CN_6737_A1615G_T{}T_gyrA'.format(i),['synonymous',['FLQ']]])

		for i in range(0,88):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrA'.format(i),[False, [False]]])

		for i in range(96,150):
			tests.append(['SNP_CN_6737_A1615G_T{}A_gyrA'.format(i),[False, [False]]])

		for test,expected in tests:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify reported as {}'.format(test, result))

	def test_classifier_COM_RIF(self):
		"""Test commercial classifier to see if it can detect positive and negative for RIF commercial mutations"""
		
		tests = []

		for i in range(424,454):
			tests.append(['SNP_CN_761141_C1335G_H{}Q_rpoB'.format(i),['asynonymous',['RIF']]])

		for i in range(424,454):
			tests.append(['SNP_CN_761141_C1335G_H{}H_rpoB'.format(i),['synonymous',['RIF']]])

		for i in range(0,424):
			tests.append(['SNP_CN_6737_A1615G_T{}A_rpoB'.format(i),[False, [False]]])

		for i in range(455,1000):
			tests.append(['SNP_CN_6737_A1615G_T{}A_rpoB'.format(i),[False, [False]]])

		for test,expected in tests:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify reported as {}'.format(test, result))

	def test_classifier_COM_INH(self):
		"""Test commercial classifier to see if it can detect positive and negative for INH commercial mutations"""
		
		tests = []
		tests.append(['SNP_CN_761141_C1335G_S315T_katG',['asynonymous',['INH']]])
		tests.append(['SNP_CN_761141_C1335G_S315S_katG',['synonymous',['INH']]])
		tests.append(['SNP_CN_761141_C1335G_S314S_katG',[False, [False]]])
		tests.append(['SNP_CN_761141_C1335G_S316S_katG',[False, [False]]])
		tests.append(['SNP_CN_761141_C1335G_S400S_katG',[False, [False]]])
		tests.append(['SNP_P_1673432_T8C_promoter-fabG1-inhA',['asynonymous',['INH']]])
		tests.append(['SNP_P_1673432_T8T_promoter-fabG1-inhA',['synonymous',['INH']]])
		tests.append(['SNP_P_1673432_T15C_promoter-fabG1-inhA',['asynonymous',['INH']]])
		tests.append(['SNP_P_1673432_T15T_promoter-fabG1-inhA',['synonymous',['INH']]])
		tests.append(['SNP_P_1673432_T16C_promoter-fabG1-inhA',['asynonymous',['INH']]])
		tests.append(['SNP_P_1673432_T16T_promoter-fabG1-inhA',['synonymous',['INH']]])
		tests.append(['SNP_P_1673432_T17C_promoter-fabG1-inhA',[False, [False]]])
		tests.append(['SNP_P_1673432_T13T_promoter-fabG1-inhA',[False, [False]]])

		for test,expected in tests:
			result = classify_COM_mutation(test)
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify as {} reported as {}'.format(test, expected, result))

	def test_classifier_WGS_INH(self):
		"""Test WGS classifier to see if it can detect AJRCCM and table-10 mutations accurately"""
		
		tests = []
		tests.append(['SNP_CN_761141_C1335G_S315T_katG',[True,'INH']])

		for test,expected in tests:
			result = commercial_WGS_tester_full.check_snp(break_down_mutation(test))
			self.assertTrue((result[0] == expected[0]) and (result[1] == expected[1]), msg = '{} failed to classify as {} reported as {}'.format(test, expected, result))

	def test_exclude_lineage_from_WGS(self):
		"""Test to see if WGS does not include any lineage snps"""
		WGS = pd.read_csv('WGS_aggregated_test_results',sep='\t')
		mutations = list(WGS['mutation'].values)
		lineage_snps_processed = ['_'.join(i.rstrip().split('_')[-2:]) for i in open('lineage_snp','r').readlines()]
		lineage_snps_unprocessed = [i.rstrip() for i in open('lineage_snp','r').readlines()]
		for mutation in mutations:
			for lineage_snp in lineage_snps_processed:
				self.assertTrue(not (lineage_snp in mutation), msg = '{} lineage snp found in {}'.format(lineage_snp, mutation))

		for mutation in mutations:
			for lineage_snp in lineage_snps_unprocessed:
				the_same = break_down_mutation(lineage_snp).compare_variant_name_location_AAchange(break_down_mutation(mutation))
				self.assertTrue(not the_same, msg = '{} lineage snp found in {}'.format(lineage_snp, mutation))


		#make sure we added those extra snps
		self.assertTrue('SNP_CN_2518919_G805A_G269S_kasA' in lineage_snps_unprocessed, msg = 'SNP_CN_2518919_G805A_G269S_kasA not found in lineage list')
		self.assertTrue('SNP_CN_2726323_C131G_P44R_ahpC' in lineage_snps_unprocessed, msg = 'SNP_CN_2726323_C131G_P44R_ahpC not found in lineage list')
		self.assertTrue('SNP_CN_2519048_G934A_G312S_kasA' in lineage_snps_unprocessed, msg = 'SNP_CN_2519048_G934A_G312S_kasA not found in lineage list')

	def test_lineage_included_commercial(self):
		"""Rest to make sure lineage snps are included in commercial results """
		commercial = pd.read_csv('commercial_aggregated_test_results',sep='\t')
		mutations = list(commercial['mutation'].values)
		lineage_snps_processed = ['_'.join(i.rstrip().split('_')[-2:]) for i in open('lineage_snp','r').readlines()]
		lineage_snps_unprocessed = [i.rstrip() for i in open('lineage_snp','r').readlines()]
		found = False
		for mutation in mutations:
			for lineage_snp in lineage_snps_processed:
				if(lineage_snp in mutation):
					found = True
			for lineage_snp in lineage_snps_unprocessed:
				if(break_down_mutation(lineage_snp).compare_variant_name_location_AAchange(break_down_mutation(mutation))):
					found = True

		self.assertTrue(not found, msg='could not find lineage snp in commercial plz check to make sure accidential exclusion did not occur')

	def test_no_synonymous_WGS(self):
		"""Test to make sure WGS test has no synonymous mutations"""
		WGS = pd.read_csv('WGS_aggregated_test_results',sep='\t')
		mutations = [break_down_mutation(i) for i in list(WGS['mutation'].values)]
		for mutation in mutations:
			type_change = mutation.AA_change
			synonymous = True
			if(len(type_change) > 2):
				#return sysnonymous since something that large probs wont be synonymous
				synonymous = False
			#If its just an insertion or something again not gonna be synonymous
			elif(len(type_change) == 1):
				synonymous = False
			elif(type_change[0] != type_change[1]):
				synonymous = False

			self.assertTrue(synonymous == False, msg='FOUND {} MUTATION WHICH IS synonymous!!!'.format(mutation))

		# self.assertTrue(WGS['synonymous'].sum() == 0, msg='synonymous mutations in WGS results')

	def test_synonymous_commercial(self):
		"""Test to make sure commercial test has synonymous mutations"""
		commercial = pd.read_csv('commercial_aggregated_test_results',sep='\t')
		mutations = [break_down_mutation(i) for i in list(commercial['mutation'].values)]
		found_synonymous= False
		for mutation in mutations:
			type_change = mutation.AA_change
			synonymous = True
			if(len(type_change) > 2):
				#return sysnonymous since something that large probs wont be synonymous
				synonymous = False
			#If its just an insertion or something again not gonna be synonymous
			elif(len(type_change) == 1):
				synonymous = False
			elif(type_change[0] != type_change[1]):
				synonymous = False
			if(not synonymous):
				found_synonymous = True

		self.assertTrue(found_synonymous == True, msg='Woops! did not find synonymous mutations in commercial make sure you did not filter them out'.format(mutation))

	def test_commercial_reclassification(self):
		"""Test to make sure all FLQ/SLI resistant ARE NOT ACCIDENTLY reclassified depending if have both INH and RIF resistance mutations"""
		WGS = pd.read_csv('commercial_raw_test_results',sep='\t')
		strains = list(set(WGS['strain'].values))
		
		for strain in strains:
			drugs_present = list(set(WGS[WGS['strain'] == strain]['drug'].values))

			#Test to make sure if no INH and no RIF that we made all of them into "IGNORE" extra annotation
			if(not('INH' in drugs_present) and not('RIF' in drugs_present)):
				if('FLQ' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'FLQ')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to FLQ but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('SLIS' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'SLIS')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to SLIS but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('EMB' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'EMB')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to EMB but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('PZA' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'PZA')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to PZA but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('STR' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'STR')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to STR but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
			else:
				#Test to make sure we did not accidently classify something as IGNORE that should not have been
				if('FLQ' in drugs_present):
					values = WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'FLQ')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to FLQ but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('SLIS' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'SLIS')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to SLIS but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('EMB' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'EMB')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to EMB but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('PZA' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'PZA')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to PZA but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('STR' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'STR')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to STR but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
		#Some more checking -- check how many should be changed vs how many were changed


	def test_WGS_reclassify_reclassification(self):
		"""Test to make sure all FLQ/SLI resistant classified isolates in WGS WGS_all tests are appropriately classified as so i.e. have both INH and RIF resistance mutations"""
		WGS = commercial_WGS_tester_full.reclassify_raw_WGS()
		strains = list(set(WGS['strain'].values))

		for strain in strains:
			drugs_present = list(set(WGS[WGS['strain'] == strain]['drug'].values))
			if(not('INH' in drugs_present) and not('RIF' in drugs_present)):
				if('FLQ' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'FLQ')]['extra_annotation'].values:
						self.assertTrue(i == 'IGNORE', msg='{} is marked resistant to FLQ but does not have INH and RIF resistance'.format(strain))
				if('SLIS' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'SLIS')]['extra_annotation'].values:
						self.assertTrue(i == 'IGNORE', msg='{} is marked resistant to SLIS but does not have INH and RIF resistance '.format(strain))
				if('EMB' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'EMB')]['extra_annotation'].values:
						self.assertTrue(i == 'IGNORE', msg='{} is marked resistant to EMB but does not have INH and RIF resistance '.format(strain))
				if('PZA' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'PZA')]['extra_annotation'].values:
						self.assertTrue(i == 'IGNORE', msg='{} is marked resistant to PZA but does not have INH and RIF resistance '.format(strain))
				if('STR' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'STR')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to STR but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
			else:
				#Test to make sure we did not accidently classify something as IGNORE that should not have been
				if('FLQ' in drugs_present):
					values = WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'FLQ')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to FLQ but does have INH or RIF resistance'.format(strain))
				if('SLIS' in drugs_present):
					values = WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'SLIS')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to SLIS but does have INH or RIF resistance'.format(strain))
				if('EMB' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'EMB')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to EMB but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('PZA' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'PZA')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to PZA but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('STR' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'STR')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to STR but does have INH or RIF resistance i.e. was reclassified!'.format(strain))

	def test_WGS_reclassification(self):
		"""Test to make sure all FLQ/SLI resistant ARE NOT ACCIDENTLY reclassified depending if have both INH and RIF resistance mutations in select WGS_test"""
		WGS = pd.read_csv('WGS_raw_test_results',sep='\t')
		strains = list(set(WGS['strain'].values))
		
		for strain in strains:
			drugs_present = list(set(WGS[WGS['strain'] == strain]['drug'].values))

			#Test to make sure if no INH and no RIF that we made all of them into "IGNORE" extra annotation
			if(not('INH' in drugs_present) and not('RIF' in drugs_present)):
				if('FLQ' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'FLQ')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to FLQ but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('SLIS' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'SLIS')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to SLIS but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('EMB' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'EMB')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to EMB but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('PZA' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'PZA')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to PZA but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
				if('STR' in drugs_present):
					for i in WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'STR')]['extra_annotation'].values:
						self.assertTrue(i != 'IGNORE', msg='{} is marked resistant to STR but does not have INH and RIF resistance i.e. it has been reclassified!'.format(strain))
			else:
				#Test to make sure we did not accidently classify something as IGNORE that should not have been
				if('FLQ' in drugs_present):
					values = WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'FLQ')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to FLQ but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('SLIS' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'SLIS')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to SLIS but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('EMB' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'EMB')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to EMB but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('PZA' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'PZA')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to PZA but does have INH or RIF resistance i.e. was reclassified!'.format(strain))
				if('STR' in drugs_present):
					values =  WGS[(WGS['strain'] == strain)&(WGS['drug'] == 'STR')]['extra_annotation'].values 
					self.assertTrue(not('IGNORE' in values), msg='{} is marked not resistant to STR but does have INH or RIF resistance i.e. was reclassified!'.format(strain))

if __name__ == '__main__':
	unittest.main()
