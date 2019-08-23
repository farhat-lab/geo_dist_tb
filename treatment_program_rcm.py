import pandas as pd 
import argparse
import os
import re
from treatment_program_rcm_helper import *



class commercial_WGS_tester():
	"""
	Class that performs commercial and WGS test on isolates located in vcf directory


	"""

	def __init__(self, vcf_directory, strain_info, results_modified_unknown):
		self.vcf_directory = vcf_directory
		self.strain_info = self.setup(strain_info)

	def __init__(self):
		self.vcf_directory = None

	def setup(self, strain_info_file):
		"""DOING INITIAL SET UP OF FILES"""
		strain_info = pd.read_csv(strain_info_file,sep='\t')

		#Added filtering to only use strains we had data on at time of analysis
		strains_we_have = [i.rstrip() for i in open('final_all_strains','r').readlines()]
		strain_info  = strain_info[strain_info['strain'].isin(strains_we_have)]
		strain_info = annotate(strain_info)

		results_modified_unknown = pd.read_csv('results_modified_unknown',sep='\t')
		strain_info = strain_info[strain_info['strain'].isin(results_modified_unknown['strain'].values)]

		return strain_info

	def break_down_mutation(self, mutation):
		"""Convert line of var file to variant class instantiation """
		if('\t' in mutation):
			mutation = mutation.split('\t')[5]
		if('&' in mutation):
			#If there is two mutations in the same position and gene we just analyze the first mutation since it has the same position and gene
			mutation = mutation.split('&')[0]
		line = mutation.replace('\n','')
		data = line.split('_')
		if(data[0] != 'SNP' and data[0] != 'INS' and data[0] != 'DEL' and data[0] != 'LSP'):
			drug = data[0]
			type_change_info = data[1:]
		else:
			drug = ''
			type_change_info = data
		# print(data)	

		if(type_change_info[0] == 'SNP' and type_change_info[1] in ['I','P','N']):
			#No real AA change info in this one with ones like SNP_I_2713795_C329T_inter-Rv2415c-eis
			gene_name, codonAA = type_change_info[4], type_change_info[3]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
			if('inter-Rv2451c-eis' in gene_name):
				#since original position relative to Rv2451c and not eis
				print("YAHOOO")
				codon_position = (type_change_info[2] - 2715332)*(-1)
				print(codon_position)
		elif(type_change_info[0] == 'LSP' and type_change_info[1] in ['CN','CS']):
			gene_name, codonAA = type_change_info[5], type_change_info[4]
			if('-' in codonAA):
				type_change, codon_position = re.sub('\d+[-]+\d+', '-', codonAA), re.findall('\d+[-]+\d+', codonAA)[0]
			else:
				type_change, codon_position = re.sub('\d+', '-', codonAA), re.findall('\d+', codonAA)[0]
		elif(type_change_info[0] == 'LSP' and type_change_info[1] in ['I','CZ']):
			gene_name, codonNT = type_change_info[4], type_change_info[3]
			type_change, codon_position = re.sub('\d+[-]+\d+','-', codonNT), re.findall('\d+[-]+\d+', codonNT)[0]
		elif(type_change_info[0] == 'SNP' ):
			gene_name, codonAA = type_change_info[5], type_change_info[4]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
		elif((type_change_info[0] == 'DEL' or type_change_info[0] == 'INS') and type_change_info[1] in ['I','P','NF','N','NI','NZ','ND']):
			gene_name, deletion = type_change_info[4], type_change_info[3]
			codon_position, type_change = re.findall('\d+', deletion)[0], re.findall('[AGCT]+', deletion)[0]
		elif(type_change_info[0] == 'DEL' or type_change_info[0] == 'INS'):
			gene_name, deletion = type_change_info[5], type_change_info[3]
			codon_position, type_change = re.findall('\d+', deletion)[0], re.findall('[AGCT]+', deletion)[0]
		else:
			raise Exception('Unknown mutation format {}'.format(type_change_info))

		return Variant(gene_name, codon_position, type_change, drug=drug)  

	# ###################SETTING UP COMMERCIAL/WGS TEST###########################

	def check_variant_commercial(self,variant):
		gene_name, codon_position, type_change = variant.gene_name, variant.codon_location, variant.AA_change
		if('-' not in codon_position):
			codon_position = int(codon_position)
		#Note negative signs are just positive
		return_variable = False
		if(gene_name == 'katG' and codon_position == 315):
			drug = 'ISONIAZID'
			return_variable = True
		elif('inhA' in gene_name and codon_position in [15, 16,8]):
			drug = 'ISONIAZID'
			return_variable = True
		elif(gene_name == 'rpoB' and codon_position in list(range(424,453))):
			drug = 'RIF'
			return_variable = True
		elif(gene_name == 'rrs' and codon_position in [1401,1402]):
			drug = 'SLIS'
			return_variable = True
		elif(gene_name != 'eis' and 'eis' in gene_name and  codon_position in [10,11,12,13,14,37]):
			drug = 'SLIS'
			return_variable = True
		#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
		elif(gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
			drug = 'FQ'
			return_variable = True
		#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
		elif(gene_name == 'gyrB' and codon_position in list(range(500,541))):
			drug = 'FQ'
			return_variable = True
		#So here if its a LSP 
		if(len(type_change) > 2 and return_variable):
			#return sysnonymous since something that large probs wont be synonymous
			return 'asynonymous', drug
		#If its just an insertion or something again not gonna be synonymous
		if(len(type_change) == 1 and return_variable):
			return 'asynonymous', drug
		if(return_variable and type_change[0] == type_change[1]):
			print("WOHOOOOO!!!!")
			return 'synonymous', drug
		elif(return_variable and type_change[0] != type_change[1]):
			return 'asynonymous', drug
		else:
			return False, False

	def check_RIF(self,gene):
		return 'rpoB' == gene

	def check_SLIS(self,gene):
		#drugs in SLIS 'KANAMYCIN' 'AMIKACIN' 'CAPREOMYCIN'
		#KAN genes: rrs, tlyA
		#AMK genes: rrs
		#CAPREOMYCIN genes: rrs

		return ('rrs' == gene ) or ('tlyA' == gene) or ('inter-eis-Rv2417c' == gene)

	def check_FQ(self,gene):
		#drugs in FQ: 'MOXIFLOXACIN 'CIPROFLOXACIN' 'OFLOXACIN'
		#MOX genes: gyrB, gyrA
		#CIPRO genes: gyrB, gyrA
		#OFLOX genes: gyrB, gyrA

		return ('gyrB' == gene) or ('gyrA' == gene)

	def check_INH(self,gene):
		condition_one = 'inhA' == gene
		condition_two = 'iniB' == gene
		condition_three = 'embB' == gene
		condition_four = 'promoter-fabG1-inhA' == gene
		condition_five = 'ahpC' == gene
		condition_six = 'promoter-ahpC' == gene
		condition_seven = 'kasA' == gene
		condition_eight = 'katG' == gene 

		return condition_one or condition_two or condition_three or condition_four or condition_five or condition_six or condition_seven or condition_eight

	def check_variant_WGS(self,variant):
		gene_name, codon_position, type_change = variant.gene_name, variant.codon_location, variant.AA_change
		if('-' not in codon_position):
			codon_position = int(codon_position)
		#Note negative signs are just positive
		return_variable = False
		if(self.check_INH(gene_name)):
			drug = 'ISONIAZID'
			return_variable = True
		elif(self.check_RIF(gene_name)):
			drug = 'RIF'
			return_variable = True
		elif(self.check_SLIS(gene_name)):
			drug = 'SLIS'
			return_variable = True
		elif(self.check_FQ(gene_name)):
			drug = 'FQ'
			return_variable = True

		#So here if its a LSP 
		if(len(type_change) > 2 and return_variable):
			#return sysnonymous since something that large probs wont be synonymous
			return 'asynonymous', drug
		#If its just an insertion or something again not gonna be synonymous
		if(len(type_change) == 1 and return_variable):
			return 'asynonymous', drug
		#ONLY CONSIDER ASYNONYMOUS
		if(return_variable and type_change[0] != type_change[1]):
			return 'asynonymous', drug
		else:
			return False, False

	def perform_analysis(self,combined, check, raw_result_file_name):

		count = 0
		result = pd.DataFrame(columns=['strain','drug','mutation','resistant','susceptible','synonymous','asynonymous'])

		for strain in list(combined['strain'].values):
			count += 1
			print("{}\t{}".format(count, len(list(combined['strain'].values))))
			#Only look at strains with data
			if(combined[combined['strain'] == strain]['NO_DATA'].item != 1):
				vcf = open(self.vcf_directory+'/'+strain+'.var','r')
				for line in vcf.readlines()[1:]:
					mutation = [i for i in line.split('\t') if 'SNP' in i or 'INS' in i or 'DEL' in i or 'LSP' in i][0]

					commercial, drug = check(self.break_down_mutation(line))
					if(commercial):
						#Make sure the strain has data for that drug
						if(combined[combined['strain'] == strain][drug].item != -1):
							if(combined[combined['strain'] == strain][drug].item != 1):
								if(commercial == 'synonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':0, 'susceptible':1, 'synonymous':1, 'asynonymous':0}
								elif(commercial == 'asynonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':0, 'susceptible':1, 'synonymous':0, 'asynonymous':1}
							else:
								if(commercial == 'synonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':1, 'susceptible':0, 'synonymous':1, 'asynonymous':0}
								elif(commercial == 'asynonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':1, 'susceptible':0, 'synonymous':0, 'asynonymous':1}

							result = result.append(to_append, ignore_index=True)
							print(to_append)

		result.to_csv(raw_result_file_name,sep='\t')
		return result

	def post_processing(self,result, name):
		lineage_snps = [i for i in open('lineage_snp','r').readlines()]
		print(result)
		print(result.columns)
		result['count'] = 1
		result = result[~result['mutation'].isin(lineage_snps)]
		result = result.groupby(['drug','mutation']).sum().reset_index().sort_values('drug',ascending=False).copy()
		result.to_csv(name, sep='\t')
		return result

	def perform_commercial_test(self):
		return self.post_processing(self.perform_analysis(self.strain_info, lambda variant: self.check_variant_commercial(variant), 'commercial_raw_test_results'), 'commercial_aggregated_test_results')

	def perform_WGS_test(self):
		return self.post_processing(self.perform_analysis(self.strain_info, lambda variant: self.check_variant_WGS(variant), 'WGS_raw_test_results'), 'WGS_aggregated_test_results')


def main():
	tester = commercial_WGS_tester('/home/lf61/lf61/mic_assemblies/46-annotate-vcfs-yasha/flatann2','strain_info.tsv','results_modified_unknown')
	tester.perform_commercial_test()
	tester.perform_WGS_test()


if __name__ == "__main__":
	main()







