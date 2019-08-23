import unittest
from treatment_program_rcm_helper import *
from treatment_program_rcm import *


break_down_mutation = lambda gene: commercial_WGS_tester().break_down_mutation(gene)

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
	# 	"""Test to make sure all isolates'"'


if __name__ == '__main__':
	unittest.main()