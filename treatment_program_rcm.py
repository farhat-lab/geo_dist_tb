import pandas as pd 
import argparse
import os
import re
from treatment_program_rcm_helper import *



class commercial_WGS_tester():
	"""
	Class that performs commercial and WGS test on isolates located in vcf directory


	"""

	def __init__(self, vcf_directory, strain_info, results_modified_unknown,lineage_snp, snp, ignore=False):
		if(not ignore):
			self.vcf_directory = vcf_directory
			self.strain_info = self.setup(strain_info)
			self.check_lineage = self.generate_lineage_snp_checker(lineage_snp)
			self.check_snp = self.generate_snp_checker(snp)


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
			gene_name, codonNT = type_change_info[4], type_change_info[3]
			if('-' in codonNT):
				type_change, codon_position = re.sub('\d+[-]+\d+','-', codonNT), re.findall('\d+[-]+\d+', codonNT)[0]
			else:
				type_change,codon_position = codonNT[0]+codonNT[len(codonNT)-1], codonNT[1:len(codonNT)-1]
			if('inter-eis-Rv2417c' in gene_name):
				#since original position relative to Rv2451c and not eis
				codon_position = str((int(type_change_info[2]) - 2715332))
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
			if('-' in codonAA):
				type_change, codon_position = re.sub('\d+[-]+\d+', '-', codonAA), re.findall('\d+[-]+\d+', codonAA)[0]
			else:
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
				drug = 'INH'
				return_variable = True
			elif(gene_name == 'promoter-fabG1-inhA' and codon_position in [15, 16,8]):
				drug = 'INH'
				return_variable = True
			elif(gene_name == 'rpoB' and codon_position in list(range(424,454))):
				drug = 'RIF'
				return_variable = True
			elif(gene_name == 'rrs' and codon_position in [1401,1402]):
				drug = 'SLIS'
				return_variable = True
			elif(gene_name == 'inter-eis-Rv2417c' and  codon_position in [10,11,12,13,14,37]):
				drug = 'SLIS'
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
			elif(gene_name == 'gyrA' and codon_position in [89,90,91,92,93,94]):
				drug = 'FLQ'
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
			elif(gene_name == 'gyrB' and codon_position in list(range(500,542))):
				drug = 'FLQ'
				return_variable = True
		else:
			#LSP deals with ranges so we see if the range intersects
			codon_position = codon_position.replace('--','-')
			start_codon_position = int(codon_position.split('-')[0])
			end_codon_position = int(codon_position.split('-')[1])
			#Note negative signs are just positive
			return_variable = False
			if(gene_name == 'katG' and (start_codon_position <= 315 <= end_codon_position)):
				drug = 'INH'
				return_variable = True
			elif(gene_name == 'promoter-fabG1-inhA' and True in [ start_codon_position <= i <= end_codon_position for i in [15, 16,8]]):
				drug = 'INH'
				return_variable = True
			elif(gene_name == 'rpoB' and True in [ start_codon_position <= i <= end_codon_position for i in list(range(424,454))]):
				drug = 'RIF'
				return_variable = True
			elif(gene_name == 'rrs' and  True in [ start_codon_position <= i <= end_codon_position for i in [1401, 1402]]):
				drug = 'SLIS'
				return_variable = True
			elif(gene_name == 'inter-eis-Rv2417c' and  True in [ start_codon_position <= i <= end_codon_position for i in [10,11,12,13,14,37]]):
				drug = 'SLIS'
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
			elif(gene_name == 'gyrA' and True in [ start_codon_position <= i <= end_codon_position for i in [89,90,91,92,93,94]]):
				drug = 'FLQ'
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
			elif(gene_name == 'gyrB' and True in [ start_codon_position <= i <= end_codon_position for i in list(range(500,542))]):
				drug = 'FLQ'
				return_variable = True


		#So here if its a LSP 
		if(len(type_change) > 2 and return_variable):
			#return sysnonymous since something that large probs wont be synonymous
			return 'asynonymous', drug
		#If its just an insertion or something again not gonna be synonymous
		if(len(type_change) == 1 and return_variable):
			return 'asynonymous', drug
		if(return_variable and type_change[0] == type_change[1]):
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
		condition_nine = 'promoter-embA-embB' == gene

		return condition_one or condition_two or condition_three or condition_four or condition_five or condition_six or condition_seven or condition_eight or condition_nine

	def check_variant_WGS(self,variant):
		gene_name, codon_position, type_change = variant.gene_name, variant.codon_location, variant.AA_change
		if('-' not in codon_position):
			codon_position = int(codon_position)
		#Note negative signs are just positive
		return_variable = False
		if(self.check_INH(gene_name)):
			drug = 'INH'
			return_variable = True
		elif(self.check_RIF(gene_name)):
			drug = 'RIF'
			return_variable = True
		elif(self.check_SLIS(gene_name)):
			drug = 'SLIS'
			return_variable = True
		elif(self.check_FQ(gene_name)):
			drug = 'FLQ'
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

	def perform_analysis(self,combined, check, raw_result_file_name, exclude_lineage, include_only_snp):
		count = 0
		result = pd.DataFrame(columns=['strain','drug','mutation','resistant','susceptible','synonymous','asynonymous', 'extra_annotation'])

		for strain in list(combined['strain'].values):
			extra_annotation = 'None'
			count += 1
			print("{}\t{}".format(count, len(list(combined['strain'].values))))
			#Only look at strains with data
			if(combined[combined['strain'] == strain]['NO_DATA'].item() != 1):
				vcf = open(self.vcf_directory+'/'+strain+'.var','r')
				for line in vcf.readlines()[1:]:
					mutation = [i for i in line.split('\t') if 'SNP' in i or 'INS' in i or 'DEL' in i or 'LSP' in i][0]
					broken_down_mutation = self.break_down_mutation(line)
					commercial, drug = check(broken_down_mutation)

					#Logic to exclude lineage mutations if needed to be detected
					if(exclude_lineage):
						if(self.check_lineage(broken_down_mutation)):
							commercial = False

					#Logic to include only 
					if(commercial and include_only_snp):
						if(self.check_snp(broken_down_mutation, drug)):
							extra_annotation = 'AJRCCM/Table-10-snp'
						else:
							extra_annotation = 'Table-10-snp'
					elif(include_only_snp):
						result, drug = self.check_snp(broken_down_mutation)
						if(result):
							commercial = True
							extra_annotation = 'AJRCCM'  
					
					if(commercial):
						#Make sure the strain has data for that drug
						if(combined[combined['strain'] == strain][drug].item() != -1):
							if(combined[combined['strain'] == strain][drug].item() != 1):
								if(commercial == 'synonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':0, 'susceptible':1, 'synonymous':1, 'asynonymous':0, 'extra_annotation':extra_annotation}
								elif(commercial == 'asynonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':0, 'susceptible':1, 'synonymous':0, 'asynonymous':1, 'extra_annotation':extra_annotation}
							else:
								if(commercial == 'synonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':1, 'susceptible':0, 'synonymous':1, 'asynonymous':0, 'extra_annotation':extra_annotation}
								elif(commercial == 'asynonymous'):
									to_append = {'strain':strain, 'drug':drug,'mutation':line.split('\t')[5], 'resistant':1, 'susceptible':0, 'synonymous':0, 'asynonymous':1, 'extra_annotation':extra_annotation}

							result = result.append(to_append, ignore_index=True)
							print(to_append)

		result.to_csv(raw_result_file_name,sep='\t')
		return result


	def generate_lineage_snp_checker(self, lineage_snp_file_location):
		#Generate function to check if a snp is a lineage snp
		lineage_snps_processed = [self.break_down_mutation(i) for i in open(lineage_snp_file_location,'r').readlines()]

		def check_if_lineage_snp(mutation):
			for lineage_snp in lineage_snps_processed:
				if(lineage_snp.compare_variant_name_location_AAchange(mutation)):
					print("{}\t{}\tFound a LINEAGE SNP".format(lineage_snp, mutation))
					return True 
			return False

		return check_if_lineage_snp


	def generate_snp_checker(self, snp_file):
		#Generate function to check if a snp is from list in AJRCCM paper
		drug_to_snp = {'INH':[],'RIF':[],'SLIS':[],'FLQ':[]}
		snp_to_drug = {}
		snps = [i.split('\t') for i in open(snp_file,'r').readlines()]
		for drug, snp in snps:
			if(drug == 'KAN' or drug == 'AMK' or drug == 'CAP'):
				drug = 'SLIS'
			elif(drug == 'LEVO' or drug == 'MOXI' or drug == 'OFLX'):
				drug = 'FLQ'
			elif(drug != 'INH' or drug != 'RIF'):
				drug = None
			if(drug):
				drug_to_snp[drug].append(self.break_down_mutation(snp.rstrip()))
				snp_to_drug[str(self.break_down_mutation(snp.rstrip()))] = drug

		def check_if_snp(mutation, drug=None):
			if(drug):	
				for snp in drug_to_snp[drug]:
					if(snp.compare_variant_name_location(mutation)):
						return True
				return False
			else:
				for snp in drug_to_snp.values():
					drug = snp_to_drug[str(mutation)]
					if(snp.compare_variant_name_location(mutation)):
						return True, drug
				return False, False

		return check_if_snp

	def post_processing(self,result, name):

		result['count'] = 1
		result = result.groupby(['drug','mutation']).sum().reset_index().sort_values('drug',ascending=False).copy()
		result.to_csv(name, sep='\t')

		return result


	def perform_commercial_test(self):
		return self.post_processing(self.perform_analysis(self.strain_info, lambda variant: self.check_variant_commercial(variant), 'commercial_raw_test_results', 0, 0), 'commercial_aggregated_test_results')

	def perform_WGS_test(self):
		return self.post_processing(self.perform_analysis(self.strain_info, lambda variant: self.check_variant_WGS(variant), 'WGS_raw_test_results', 1, 1), 'WGS_aggregated_test_results')

	def calculate_statistics_commercial(self):
		df = pd.read_csv('commercial_raw_test_results',sep='\t')
		for drug in ['ISONIAZID','RIF','SLIS','FQ']:
			number_resistant_predicted = df[(df['drug'] == drug) & (df['resistant'] == 1)]['strain'].nunique()
			number_susceptible_predicted = df[(df['drug'] == drug) & (df['susceptible'] == 1)]['strain'].nunique()

			
			number_resistant = self.strain_info[self.strain_info[drug] == 1]['strain'].count()
			number_susceptible = self.strain_info[self.strain_info[drug] == 0]['strain'].count()

			print("{} NUMBER RESISTANT NOT ABLE TO BE PREDICTED {}/{} {}".format(drug, number_resistant - number_resistant_predicted, number_resistant, 1-(number_resistant_predicted/number_resistant)))
			print("{} SUSCEPTIBLE PREDICTED {}/{} {}".format(drug, number_susceptible_predicted, number_susceptible, number_susceptible_predicted/number_susceptible))

	def calculate_statistics_WGS(self):
		df = pd.read_csv('WGS_raw_test_results',sep='\t')
		for drug in ['ISONIAZID','RIF','SLIS','FQ']:
			number_resistant_predicted = df[(df['drug'] == drug) & (df['resistant'] == 1)]['strain'].nunique()
			number_susceptible_predicted = df[(df['drug'] == drug) & (df['susceptible'] == 1)]['strain'].nunique()

			
			number_resistant = self.strain_info[self.strain_info[drug] == 1]['strain'].count()
			number_susceptible = self.strain_info[self.strain_info[drug] == 0]['strain'].count()

			print("{} NUMBER RESISTANT NOT ABLE TO BE PREDICTED {}/{} {}".format(drug, number_resistant - number_resistant_predicted, number_resistant, 1-(number_resistant_predicted/number_resistant)))
			print("{} SUSCEPTIBLE PREDICTED {}/{} {}".format(drug, number_susceptible_predicted, number_susceptible, number_susceptible_predicted/number_susceptible))

def main():
	tester = commercial_WGS_tester('/home/lf61/lf61/mic_assemblies/46-annotate-vcfs-yasha/flatann2','strain_info.tsv','results_modified_unknown', 'lineage_snp', 'snps')
	tester.perform_commercial_test()
	tester.perform_WGS_test()
	# print("COMMERCIAL STATS")
	# tester.calculate_statistics_commercial()
	# print("WGS STATS \n\n\n\n")
	# tester.calculate_statistics_WGS()


if __name__ == "__main__":
	main()







