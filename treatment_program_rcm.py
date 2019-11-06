import pandas as pd 
import argparse
import os
import re
from treatment_program_rcm_helper import *
import scipy.stats as stats


class commercial_WGS_tester():
	"""
	Class that performs commercial and WGS test on isolates located in vcf directory


	"""

	def __init__(self, vcf_directory, strain_info, results_modified_unknown,lineage_snp, snp, cryptic_snp, ignore=False):
		if(not ignore):
			self.vcf_directory = vcf_directory
			self.strain_info = self.setup(strain_info)
			self.check_lineage = self.generate_lineage_snp_checker(lineage_snp)
			self.check_snp = self.generate_snp_checker(snp)
			self.check_snp_cryptic = self.generate_cryptic_snp_checker(cryptic_snp)


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
		elif((type_change_info[0] == 'DEL' or type_change_info[0] == 'INS') and type_change_info[1] in ['I','F','P','NF','N','NI','NZ','ND']):
			gene_name, deletion = type_change_info[4], type_change_info[3]
			if('i' not in deletion and 'd' not in deletion):
				raise Exception('VARIANT {} does not have i or d which we assume is present'.format(line))
			codon_position, type_change = re.findall('\d+', deletion)[0], re.findall('[AGCT]+', deletion)[0]
		elif(type_change_info[0] == 'DEL' or type_change_info[0] == 'INS'):
			gene_name, deletion = type_change_info[5], type_change_info[3]
			if('i' not in deletion and 'd' not in deletion):
				raise Exception('VARIANT {} does not have i or d which we assume is present'.format(line))
			codon_position, type_change = re.findall('\d+', deletion)[0], re.findall('[AGCT]+', deletion)[0]
		else:
			raise Exception('Unknown mutation format {}'.format(type_change_info))

		return Variant(gene_name, codon_position, type_change, line, drug=drug)  

	# ###################SETTING UP COMMERCIAL/WGS TEST###########################

	def check_variant_commercial(self,variant):
		gene_name, codon_position, type_change = variant.gene_name, variant.codon_location, variant.AA_change
		drug = []
		if('-' not in codon_position):
			codon_position = int(codon_position)
			#Note negative signs are just positive
			return_variable = False
			if(gene_name == 'katG' and codon_position == 315):
				drug.append('INH')
				return_variable = True
			elif(gene_name == 'promoter-fabG1-inhA' and codon_position in [15, 16,8]):
				drug.append('INH')
				return_variable = True
			if(gene_name == 'rpoB' and codon_position in list(range(424,454))):
				drug.append('RIF')
				return_variable = True
			if(gene_name == 'rrs' and codon_position in [1401,1402]):
				drug.append('SLIS')
				return_variable = True
			elif(gene_name == 'inter-eis-Rv2417c' and  codon_position in [10,11,12,13,14,37]):
				drug.append('SLIS')
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
			if(gene_name == 'gyrA' and codon_position in [89,90,91,92,93,94]):
				drug.append('FLQ')
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
			elif(gene_name == 'gyrB' and codon_position in list(range(500,542))):
				drug.append('FLQ')
				return_variable = True
		else:
			#LSP deals with ranges so we see if the range intersects
			codon_position = codon_position.replace('--','-')
			start_codon_position = int(codon_position.split('-')[0])
			end_codon_position = int(codon_position.split('-')[1])
			#Note negative signs are just positive
			return_variable = False
			if(gene_name == 'katG' and (start_codon_position <= 315 <= end_codon_position)):
				drug.append('INH')
				return_variable = True
			elif(gene_name == 'promoter-fabG1-inhA' and True in [ start_codon_position <= i <= end_codon_position for i in [15, 16,8]]):
				drug = 'INH'
				return_variable = True
			if(gene_name == 'rpoB' and True in [ start_codon_position <= i <= end_codon_position for i in list(range(424,454))]):
				drug.append('RIF')
				return_variable = True
			if(gene_name == 'rrs' and  True in [ start_codon_position <= i <= end_codon_position for i in [1401, 1402]]):
				drug.append('SLIS')
				return_variable = True
			elif(gene_name == 'inter-eis-Rv2417c' and  True in [ start_codon_position <= i <= end_codon_position for i in [10,11,12,13,14,37]]):
				drug.append('SLIS')
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
			if(gene_name == 'gyrA' and True in [ start_codon_position <= i <= end_codon_position for i in [89,90,91,92,93,94]]):
				drug.append('FLQ')
				return_variable = True
			#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
			elif(gene_name == 'gyrB' and True in [ start_codon_position <= i <= end_codon_position for i in list(range(500,542))]):
				drug.append('FLQ')
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
			return False, [False]

	def check_RIF(self,gene):
		return 'rpoB' == gene

	def check_STR(self, gene):
		return ('gid' == gene) or ('rpsL' == gene) or ('rrs' == gene) or ('inter-rrs-rrl' == gene)

	def check_PZA(self, gene):
		return ('promoter-pncA' == gene) or ('pncA' == gene) or ('rpsA' == gene)

	def check_EMB(self, gene):
		return ('embA' == gene) or ('embB' == gene) or ('embC' == gene) or ('iniB' == gene) or ('promoter-embA-embB' == gene)

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
		condition_ten = 'fabG1' == gene

		return condition_one or condition_two or condition_three or condition_four or condition_five or condition_six or condition_seven or condition_eight or condition_nine or condition_ten

	def check_variant_WGS(self,variant):
		gene_name, codon_position, type_change = variant.gene_name, variant.codon_location, variant.AA_change
		if('-' not in codon_position):
			codon_position = int(codon_position)
		#Note negative signs are just positive
		drug = []
		return_variable = False
		if(self.check_INH(gene_name)):
			drug.append('INH')
			return_variable = True
		if(self.check_RIF(gene_name)):
			drug.append('RIF')
			return_variable = True
		if(self.check_SLIS(gene_name)):
			drug.append('SLIS')
			return_variable = True
		if(self.check_FQ(gene_name)):
			drug.append('FLQ')
			return_variable = True
		if(self.check_PZA(gene_name)):
			drug.append('PZA')
			return_variable = True
		if(self.check_EMB(gene_name)):
			drug.append('EMB')
			return_variable = True
		if(self.check_STR(gene_name)):
			drug.append('STR')
			return_variable = True
		#So here if its a LSP 
		if(len(type_change) > 2 and return_variable):
			#return sysnonymous since something that large probs wont be synonymous
			variant.synonymous = False
			return 'asynonymous', drug
		#If its just an insertion or something again not gonna be synonymous
		if(len(type_change) == 1 and return_variable):
			variant.synonymous = False
			return 'asynonymous', drug
		#ONLY CONSIDER ASYNONYMOUS
		if(return_variable and type_change[0] != type_change[1]):
			variant.synonymous = False
			return 'asynonymous', drug
		else:
			variant.synonymous = True
			return False, [False]

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
					commercial, drugs = check(broken_down_mutation)
					#Logic to exclude lineage mutations if needed to be detected
					no_lineage = True
					if(exclude_lineage):
						#This means we are in a WGS test were being asynonymous and no lineage snp matters
						if(self.check_lineage(broken_down_mutation)):
							commercial = False
							no_lineage = False
					else:
						#We are in a commercial test where synonymous v asynonymous and lineage doesn't matter
						break_down_mutation.synonymous = False
					old_commercial = commercial
					for drug in drugs:
						commercial = old_commercial
						#Logic to include only 
						if(commercial and include_only_snp):
							AJRCCM_result = self.check_snp(broken_down_mutation, drug)
							CRYPTIC_result = self.check_snp_cryptic(broken_down_mutation, drug)
							if(AJRCCM_result):
								if(CRYPTIC_result):
									extra_annotation = 'AJRCCM/CRYPTIC/Table-10-snp'
								else:
									extra_annotation = 'AJRCCM/Table-10-snp'
							elif(CRYPTIC_result):
								extra_annotation = 'CRYPTIC/Table-10-snp'
							else:
								extra_annotation = 'Table-10-snp'
						elif(include_only_snp and not broken_down_mutation.synonymous and no_lineage):
							AJRCCM_test_result, AJRCCM_drug = self.check_snp(broken_down_mutation)
							CRYPTIC_test_result, CRYPTIC_drug = self.check_snp_cryptic(broken_down_mutation)
							if(AJRCCM_test_result):
								if(CRYPTIC_test_result):
									commercial = True
									extra_annotation = 'CRYPTIC/AJRCCM'
								else:
									commercial = True
									extra_annotation = 'AJRCCM'  
							elif(CRYPTIC_test_result):
								commercial = True
								extra_annotation = 'CRYPTIC'
							if(commercial):
								print(AJRCCM_test_result)
								print(CRYPTIC_test_result)
								raise Exception("ALL WGS test failed on something that select WGS succeeded {}".format(line))
					
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
		# #RECLASSIFICATION
		# result = self.perform_reclassification(result)

		result.to_csv(raw_result_file_name,sep='\t')
		return result


	def generate_lineage_snp_checker(self, lineage_snp_file_location):
		#Generate function to check if a snp is a lineage snp
		lineage_snps_processed = [self.break_down_mutation(i) for i in open(lineage_snp_file_location,'r').readlines()]

		def check_if_lineage_snp(mutation):
			for lineage_snp in lineage_snps_processed:
				if(lineage_snp.compare_variant_name_location_AAchange(mutation)):
					return True 
			return False

		return check_if_lineage_snp


	def generate_snp_checker(self, snp_file):
		#Generate function to check if a snp is from list in AJRCCM paper
		drug_to_snp = {'INH':[],'RIF':[],'SLIS':[],'FLQ':[], 'PZA':[], 'STR':[], 'EMB':[]}
		snp_to_drug = {}
		snps = [i.split('\t') for i in open(snp_file,'r').readlines()]
		for drug, snp in snps:
			if(drug == 'KAN' or drug == 'AMK' or drug == 'CAP'):
				drug = 'SLIS'
			elif(drug == 'LEVO' or drug == 'MOXI' or drug == 'OFLX'):
				drug = 'FLQ'
			elif(drug != 'INH' and drug != 'RIF' and drug != 'PZA' and drug != 'STR' and drug != 'EMB'):
				drug = None
			if(drug):
				drug_to_snp[drug].append(self.break_down_mutation(snp.rstrip()))
				snp_to_drug[str(self.break_down_mutation(snp.rstrip()))] = drug
		def check_if_snp(mutation, drug=None):
			if(drug):	
				for snp in drug_to_snp[drug]:
					if(snp.compare_variant_name_location(mutation)):
						return True
				if(drug == 'RIF'):
					if(mutation.gene_name == 'rpoB' and mutation.is_frameshift):
						return True
				elif(drug == 'PZA'):
					if(mutation.gene_name == 'pncA' and mutation.is_frameshift):
						return True
				elif(drug == 'INH'):
					if(mutation.gene_name == 'katG' and mutation.is_frameshift):
						return True
				return False
			else:
				for snp in [item for sublist in drug_to_snp.values() for item in sublist]:
					if(snp.compare_variant_name_location(mutation)):
						drug = snp_to_drug[str(mutation)]
						return True, drug
				if(mutation.gene_name == 'rpoB' and mutation.is_frameshift):
					return True, 'RIF'
				elif(mutation.gene_name == 'pncA' and mutation.is_frameshift):
					return True, 'PZA'
				elif(mutation.gene_name == 'katG' and mutation.is_frameshift):
					return True, 'INH'
				return False, False

		return check_if_snp

	def cryptic_snp_parser(self, line):
		gene_name, change, type_change = [i.rstrip() for i in line.split('\t')]
		drug = []

		if(type_change == 'PROT'): 
			specific_change,codon_position = change[0]+change[len(change)-1], change[1:len(change)-1]
		elif(type_change == 'DNA'):
			if(len(change) > 10):
				#Then this is a LSP
				if(gene_name == 'pncA'):
					gene_name = 'promoter-pncA'
				else:
					raise Exception("gene name we do not have promoter to {}".format(line))
				
				specific_change, codon_position = re.sub('[-]+\d+', '-', change), re.findall('\d+', change)[0]
			else:
				if('-' in change):
					if(gene_name == 'eis'):
						gene_name = 'inter-eis-Rv2417c'
					elif(gene_name == 'pncA'):
						gene_name = 'promoter-pncA'
					elif(gene_name == 'embA'):
						gene_name = 'promoter-embA-embB'
					elif(gene_name == 'ahpC'):
						gene_name = 'promoter-ahpC'
					elif(gene_name == 'fabG1'):
						gene_name = 'promoter-fabG1-inhA'
					else:
						raise Exception("gene name we do not have promoter to {}".format(line))
					change = re.sub('-','', change)
				if(len(re.sub('\d+','',change)) > 2):
					specific_change,codon_position = re.sub('\d+', '-', change), re.findall('\d+',change)[0]
				else:
					specific_change,codon_position = change[0]+change[len(change)-1], change[1:len(change)-1]
		else:
			raise Exception("Do not recognize type {}".format(type_change))

		for potential_drug in drug_to_WGSregion:
			if(gene_name in drug_to_WGSregion[potential_drug]):
				drug.append(potential_drug) 
		if(drug == []):
			raise Exception("Could not find a drug for the gene {}".format(line))

		return drug, Variant(gene_name, codon_position, specific_change, line)

	def generate_cryptic_snp_checker(self, snp_file):
		#Generate function to check if a snp is from list in CRYPTIC snp list
		drug_to_snp = {'INH':[],'RIF':[],'SLIS':[],'FLQ':[], 'PZA':[], 'STR':[], 'EMB':[]}
		snp_to_drug = {}
		for line in open(snp_file,'r').readlines():
			drugs, mutation = self.cryptic_snp_parser(line)
			for drug in drugs:
				drug_to_snp[drug].append(mutation)
				snp_to_drug[str(mutation)] = drug

		def check_if_snp(mutation, drug=None):
			if(drug):	
				for snp in drug_to_snp[drug]:
					if(snp.compare_variant_name_location(mutation)):
						return True
				if(drug == 'RIF'):
					if(mutation.gene_name == 'rpoB' and mutation.is_frameshift):
						return True
				elif(drug == 'PZA'):
					if(mutation.gene_name == 'pncA' and mutation.is_frameshift):
						return True
				elif(drug == 'INH'):
					if(mutation.gene_name == 'katG' and mutation.is_frameshift):
						return True
				return False
			else:
				for snp in [item for sublist in drug_to_snp.values() for item in sublist]:
					if(snp.compare_variant_name_location(mutation)):
						drug = snp_to_drug[str(mutation)]
						return True, drug
				if(mutation.gene_name == 'rpoB' and mutation.is_frameshift):
					return True, 'RIF'
				elif(mutation.gene_name == 'pncA' and mutation.is_frameshift):
					return True, 'PZA'
				elif(mutation.gene_name == 'katG' and mutation.is_frameshift):
					return True, 'INH'
				return False, False

		return check_if_snp


	def post_processing(self,result, name, recreate = False):
		if(recreate):
			result = pd.read_csv(recreate, sep='\t')

		result['count'] = 1
		result = result.groupby(['drug','mutation', 'extra_annotation']).sum().reset_index().sort_values('drug',ascending=False).copy()
		#del result['Unnamed: 0']
		result.to_csv(name, sep='\t')

		return result

	def perform_reclassification(self, result):
		strains = list(set(result['strain'].values))
		for strain in strains:
			if(self.SLIFLQ_reclassification_required(strain, result)):
				#We change all resistant to 0 and susceptible to 1 and type to "IGNORE" to ignore the mutation
				result.loc[(result.strain==strain) & (result['drug'].isin(['FLQ','SLIS', 'PZA','EMB'])), 'resistant'] = 0
				result.loc[(result.strain==strain) & (result['drug'].isin(['FLQ','SLIS','PZA','EMB'])), 'susceptible'] = 1
				result.loc[(result.strain==strain) & (result['drug'].isin(['FLQ','SLIS','PZA','EMB'])), 'extra_annotation'] = 'IGNORE'
				# print("{} RECLASSIFYING".format(strain))
			#else:
				# print("{} IS GOOD".format(strain))

		return result


	def SLIFLQ_reclassification_required(self, strain, result):
		drugs_present = list(set(result[result['strain'] == strain]['drug'].values))
		#was the strain marked resistance to FLQ or SLI?
		return not('INH' in drugs_present )  and not ('RIF' in drugs_present)
		# if(result[(result['strain'] == strain)&(result['drug'] == 'FLQ')]['resistant'].sum().item() or result[(result['strain'] == strain)&(result['drug'] == 'SLIS')]['resistant'].sum().item()):
		# 	return not ('INH' in drugs_present and 'RIF' in drugs_present)

	def perform_commercial_test(self):
		return self.post_processing(self.perform_analysis(self.strain_info, lambda variant: self.check_variant_commercial(variant), 'commercial_raw_test_results', 0, 0), 'commercial_aggregated_test_results')

	def perform_WGS_test(self):
		return self.post_processing(self.perform_analysis(self.strain_info, lambda variant: self.check_variant_WGS(variant), 'WGS_raw_test_results', 1, 1), 'WGS_aggregated_test_results')

	def perform_commercial_post_processing(self):
		return self.post_processing(None, 'commercial_aggregated_test_results', 'commercial_raw_test_results')

	def perform_WGS_post_processing(self):
		return self.post_processing(None, 'WGS_aggregated_test_results', 'WGS_raw_test_results')

	def reclassify_raw_result(self, name):
		result = pd.read_csv(name,sep='\t')
		#self.perform_reclassification(result).to_csv(name,sep='\t')
		return self.perform_reclassification(result)

	def reclassify_raw_commercial(self):
		self.reclassify_raw_result('commercial_raw_test_results')

	def reclassify_raw_WGS(self):
		return self.reclassify_raw_result('WGS_raw_test_results')

	# def get_revised_strain_info_count(self, modified_raw_test_results, drug):
	# 	"""After we revise phenotype of the isolates we need to count via strain info but make sure phenotype wasn't revised in our analysis"""
	# 	resistant_count = 0
	# 	susceptible_count = 0
	# 	#First we cycle through all relevant strains for the drug that have data
	# 	for strain in self.strain_info[self.strain_info[drug] != -1]['strain'].values:
	# 		#If the strain is resistant to the drug
	# 		if(self.strain_info[(self.strain_info['strain'] ==strain)][drug].item()):
	# 			#We check if the strain is still resistant to that drug i.e. it wasn't reclassified
	# 			if(modified_raw_test_results[(modified_raw_test_results['strain'] == strain)&(modified_raw_test_results['drug'] == drug)]['resistant'].sum().item()):
	# 				#If this is true we up the resistant count
	# 				resistant_count += 1
	# 			else:
	# 				#If this is not true we reclassified this as susceptible so we up the susceptible count
	# 				susceptible_count += 1

	# 	#Now we only reclassified resistant isolates, susceptible ones stayed the same so we we can just add the same as we used to do it
	# 	susceptible_count += self.strain_info[self.strain_info[drug] == 0]['strain'].count()

		# return resistant_count, susceptible_count


	def calculate_statistics_commercial(self):
		df = pd.read_csv('commercial_raw_test_results',sep='\t')
		df = df[df['extra_annotation'] != 'IGNORE']
		for drug in ['INH','RIF','SLIS','FLQ']:
			number_resistant_predicted = df[(df['drug'] == drug) & (df['resistant'] == 1)]['strain'].nunique()
			number_susceptible_predicted = df[(df['drug'] == drug) & (df['susceptible'] == 1)]['strain'].nunique()

			
			# number_resistant, number_susceptible = self.get_revised_strain_info_count(df, drug)
			number_resistant = self.strain_info[self.strain_info[drug] == 1]['strain'].count()
			number_susceptible = self.strain_info[self.strain_info[drug] == 0]['strain'].count()

			print("{} NUMBER RESISTANT PREDICTED (SENSITIVITY) {}/{} {}".format(drug, number_resistant_predicted, number_resistant, float(number_resistant_predicted)/number_resistant))
			print("{} SUSCEPTIBLE CORRECTLY PREDICTED AS NOT RESISTANT (SPECIFICITY) {}/{} {}".format(drug, number_susceptible - number_susceptible_predicted, number_susceptible, 1-(number_susceptible_predicted/float(number_susceptible))))

		print("CALCULATING SENSITIVITY/SPECIFICTY ACROSS COUNTRIES WITH TOP 5 AMOUNT DATA")
		#First determine country with top 5 amount of data
		df = pd.merge(df, self.strain_info[['strain','country']], on='strain',how='inner')
		top_five_countries = [i[0] for i in df[df['country'] != 'Not Provided'].groupby('country').country.value_counts().nlargest(5).index.tolist()]
		for country in top_five_countries:
			print("STATISTICS FOR {}".format(country))
			local_df = df[df['country'] == country].copy()
			for drug in ['INH','RIF','SLIS','FLQ']:
				number_resistant_predicted = local_df[(local_df['drug'] == drug) & (local_df['resistant'] == 1)]['strain'].nunique()
				number_susceptible_predicted = local_df[(local_df['drug'] == drug) & (local_df['susceptible'] == 1)]['strain'].nunique()

				# number_resistant, number_susceptible = self.get_revised_strain_info_count(df, drug)
				number_resistant = self.strain_info[(self.strain_info[drug] == 1) & (self.strain_info['country'] == country)]['strain'].count()
				number_susceptible = self.strain_info[(self.strain_info[drug] == 0) & (self.strain_info['country'] == country)]['strain'].count()
				
				if(number_resistant == 0):
					print("{} NUMBER RESISTANT PREDICTED (SENSITIVITY) {}/{} {}".format(drug, number_resistant_predicted, number_resistant, "0/0 NA"))
				else:
					print("{} NUMBER RESISTANT PREDICTED (SENSITIVITY) {}/{} {}".format(drug, number_resistant_predicted, number_resistant, float(number_resistant_predicted)/number_resistant))
				
				if(number_susceptible == 0):
					print("{} SUSCEPTIBLE CORRECTLY PREDICTED AS NOT RESISTANT (SPECIFICITY) {}/{} {}".format(drug, number_susceptible - number_susceptible_predicted, number_susceptible, "0/0 NA" ))
				else:
					print("{} SUSCEPTIBLE CORRECTLY PREDICTED AS NOT RESISTANT (SPECIFICITY) {}/{} {}".format(drug, number_susceptible - number_susceptible_predicted, number_susceptible, 1-(number_susceptible_predicted/float(number_susceptible))))
		
		#Calculate P-value test for Peru v South Africa Sensitivity for SLI and FLQ
		for drug, country_one, country_two in [['SLIS','Peru','South Africa'],['FLQ','Peru','South Africa']]:
			local_df = df[df['country'] == country_one].copy()
			number_resistant_predicted_one = local_df[(local_df['drug'] == drug) & (local_df['resistant'] == 1)]['strain'].nunique()
			number_resistant_one =  self.strain_info[(self.strain_info[drug] == 1) & (self.strain_info['country'] == country_one)]['strain'].count()

			local_df = df[df['country'] == country_two].copy()
			number_resistant_predicted_two = local_df[(local_df['drug'] == drug) & (local_df['resistant'] == 1)]['strain'].nunique()
			number_resistant_two =  self.strain_info[(self.strain_info[drug] == 1) & (self.strain_info['country'] == country_two)]['strain'].count()

			table = [[number_resistant_one, number_resistant_one-number_resistant_predicted_one], 
			[number_resistant_two, number_resistant_two-number_resistant_predicted_two]]
			print("TESTS FOR {}".format(drug))
			print("{}:\t{}/{} {}".format(country_one, number_resistant_predicted_one, number_resistant_one, number_resistant_predicted_one/float(number_resistant_one)))
			print("{}:\t{}/{} {}".format(country_two, number_resistant_predicted_two, number_resistant_two, number_resistant_predicted_two/float(number_resistant_two)))

			oddsratio, pvalue = stats.fisher_exact(table)
			print(pvalue)

		


	def calculate_statistics_WGS(self):
		# total_df = pd.read_csv('WGS_raw_test_results',sep='\t')
		# total_df = total_df[total_df['extra_annotation'] != 'IGNORE']
		for annotation in ['Table-10-snp','AJRCCM', 'CRYPTIC']:
			print("FOR FOLLOWING MUTATIONS IN {}".format(annotation))
			if(annotation == 'Table-10-snp'):
				total_df = self.reclassify_raw_WGS()
				total_df = total_df[total_df['extra_annotation'] != 'IGNORE']
			else:
				total_df = pd.read_csv('WGS_raw_test_results',sep='\t')

			df = total_df[total_df['extra_annotation'].str.contains(annotation)]
			for drug in ['INH','RIF','SLIS','FLQ', 'STR','PZA','EMB']:
				number_resistant_predicted = df[(df['drug'] == drug) & (df['resistant'] == 1)]['strain'].nunique()
				number_susceptible_predicted = df[(df['drug'] == drug) & (df['susceptible'] == 1)]['strain'].nunique()

				# number_resistant, number_susceptible = self.get_revised_strain_info_count(df, drug)
				number_resistant = self.strain_info[self.strain_info[drug] == 1]['strain'].count()
				number_susceptible = self.strain_info[self.strain_info[drug] == 0]['strain'].count()

				print("{} NUMBER RESISTANT PREDICTED (SENSITIVITY) {}/{} {}".format(drug, number_resistant_predicted, number_resistant, float(number_resistant_predicted)/number_resistant))
				print("{} SUSCEPTIBLE CORRECTLY PREDICTED AS NOT RESISTANT (SPECIFICITY) {}/{} {}".format(drug, number_susceptible - number_susceptible_predicted, number_susceptible, 1-(number_susceptible_predicted/float(number_susceptible))))

def main():
	#NOTE ONLY RECLASSIFICATION HAPPENS FOR ALL_WGS_TEST!
	tester = commercial_WGS_tester('/home/lf61/lf61/mic_assemblies/46-annotate-vcfs-yasha/flatann2','strain_info.tsv','results_modified_unknown', 'lineage_snp', 'snps', 'panel.final.Cryptic_no_frameshift.txt')
	#tester.perform_commercial_test()
	# tester.perform_WGS_test()
	#tester.perform_commercial_post_processing()
	# tester.perform_WGS_post_processing()
	# tester.reclassify_raw_commercial()
	# tester.reclassify_raw_WGS()
	print("COMMERCIAL STATS")
	tester.calculate_statistics_commercial()
	print("WGS STATS \n\n\n\n")
	tester.calculate_statistics_WGS()
	
	


if __name__ == "__main__":
	main()







