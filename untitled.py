import unittest
from treatment_program_rcm_helper import *
from treatment_program_rcm import *

commercial_WGS_tester_full = commercial_WGS_tester('/home/lf61/lf61/mic_assemblies/46-annotate-vcfs-yasha/flatann2','strain_info.tsv','results_modified_unknown', 'lineage_snp', 'snps','panel.final.Cryptic_no_frameshift.txt')
commercial_WGS_tester = commercial_WGS_tester('ugh','ugh','ugh','ugh','ugh','ugh',ignore=True)
break_down_mutation = lambda gene: commercial_WGS_tester.break_down_mutation(gene)
classify_COM_mutation =  lambda gene: commercial_WGS_tester.check_variant_commercial(break_down_mutation(gene))
classify_WGS_mutation = lambda gene: commercial_WGS_tester.check_variant_WGS(gene)


class Test(unittest.TestCase):



if __name__ == '__main__':
	unittest.main()
