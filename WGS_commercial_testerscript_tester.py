import unittest
from treatment_program_rcm_helper import drug_to_region

class Test(unittest.TestCase):

	def test_regions(self):
		WGS = pd.read_csv('WGS_test_results',sep='\t')
		mutation_values WGS['mutation'].values

if __name__ == '__main__':
	unittest.main()