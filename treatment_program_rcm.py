import pandas as pd 
import argparse
import os
import re
from treatment_program_rcm_helper import *

parser = argparse.ArgumentParser(description='Creates statistics about resistance detectable by commercial diagnostics vs WGS')

parser.add_argument('vcf_files', help='directory with vcf files')
parser.add_argument('strain_info_file', help='the strain_info_file with resistance data for every strain')
parser.add_argument('results_modified_unknown', help='file with output of vcfmatchunk.py')

args = parser.parse_args()
results_modified_unknown, vcf_directory = open(args.results_modified_unknown, 'r'), args.vcf_files

######DOING INITIAL SET UP OF FILES#############
strain_info = pd.read_csv(args.strain_info_file,sep='\t')

#Added filtering to only use strains we had data on at time of analysis
strains_we_have = [i.rstrip() for i in open('final_all_strains','r').readlines()]
strain_info  = strain_info[strain_info['strain'].isin(strains_we_have)]
strain_info = annotate(strain_info)

resistance_mutation = pd.read_csv('results_modified_unknown',sep='\t')
combined = pd.merge(strain_info, resistance_mutation, on='strain', how='outer')

cols = list(combined.columns)
columns_of_interest = []
for i in cols:
	for drug in drugs:
		if (drug in i and i != drug_mapping[drug]):
			columns_of_interest.append(i)

columns_of_interest.append('country')

combined = combined.dropna(subset=columns_of_interest)
del combined['Unnamed: 0']

vcf_files = [i for i in os.listdir(vcf_directory) if '.var' in i]

# ###################SETTING UP DRUG_TO_VARIANT###########################


drug_to_variants = {}

def break_down_mutation(mutation):
	if('\t' in mutation):
		mutation = mutation.split('\t')[5]
	if('&' in mutation):
		#If there is two mutations in the same position and gene we just analyze the first mutation since it has the same position and gene
		mutation = mutation.split('&')[0]
	line = mutation.replace('\n','')
	data = line.split('_')
	print(data)
	if(data[0] != 'SNP' and data[0] != 'INS' and data[0] != 'DEL' and data[0] != 'LSP'):
		drug = data[0]
		type_change_info = data[1:]
	else:
		drug = ''
		type_change_info = data
	

	if(type_change_info[0] == 'SNP' and type_change_info[1] in ['I','P','N']):
		#No real AA change info in this one with ones like SNP_I_2713795_C329T_inter-Rv2415c-eis
		gene_name, codonAA = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
	elif(type_change_info[0] == 'LSP' and type_change_info[1] in ['CN','CS']):
		gene_name, codonAA, codonNT = type_change_info[5], type_change_info[4], type_change_info[3]
		type_change, codon_position = re.sub('\d+[-]+\d+', '-', codonAA), re.findall('\d+[-]+\d+', codonNT)[0]
	elif(type_change_info[0] == 'LSP' and type_change_info[1] in ['I','CZ']):
		gene_name, codonNT = type_change_info[4], type_change_info[3]
		type_change, codon_position = re.sub('\d+[-]+\d+','-', codonNT), re.findall('\d+[-]+\d+', codonNT)[0]
	elif(type_change_info[0] == 'SNP' ):
		gene_name, codonAA, codonNT = type_change_info[5], type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonNT[1:len(codonNT)-1]
	elif((type_change_info[0] == 'DEL' or type_change_info[0] == 'INS') and type_change_info[1] in ['I','P','NF']):
		gene_name, deletion = type_change_info[4], type_change_info[3]
		codon_position, type_change = re.findall('\d+', deletion)[0], re.findall('[AGCT]+', deletion)[0]
	elif(type_change_info[0] == 'DEL' or type_change_info[0] == 'INS'):
		gene_name, deletion = type_change_info[5], type_change_info[3]
		codon_position, type_change = re.findall('\d+', deletion)[0], re.findall('[AGCT]+', deletion)[0]
	else:
		raise Exception('Unknown mutation format {}'.format(type_change_info))

	# if(type_change_info[0] == 'SNP' and type_change_info[1] == 'CN'):
	#     	gene_name, codonAA = type_change_info[5], type_change_info[4]
	# elif(type_change_info[0] == 'SNP' and type_change_info[1] == 'P'):
	#     	gene_name, codonAA = type_change_info[4], type_change_info[3]
	#     	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
	# elif(type_change_info[0] == 'SNP' and type_change_info[1] == 'N'):
	#     	gene_name, codonAA = type_change_info[4], type_change_info[3]
	#     	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
	# elif(type_change_info[0] == 'SNP' and type_change_info[1] == 'I'):
	#     	gene_name, codonAA= type_change_info[4], type_change_info[3]
	#     	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','')
	# elif(type_change_info[0] == 'INS'):
	#     	gene_name, codonAA = type_change_info[4], type_change_info[3]
	#     	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','')
	# elif(type_change_info[0] == 'DEL'):
	#     	gene_name, codonAA = type_change_info[4], type_change_info[3]
	#     	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','')
	# else:
	# 	#So this is a non-standard type, but it may still be relevant so we can just use regex
	# 	#for gene name we can always assume its in the same position
	# 	gene_name = type_change_info[4]
	# 	#However, for the codon position we use regex any numbers will be codon position
	# 	codon_position = [i for i in re.findall('\d*', type_change_info[3]) if i][0]
	# 	type_change = ''
	print("{}\t{}\t{}\t{}".format(gene_name, codon_position, type_change, drug))
	return Variant(gene_name, codon_position, type_change, drug=drug)  

###RUN CHECK MUTATION TEST i.e. see no failures
count = 0
for vcf_file in vcf_files:
	vcf = open(vcf_directory+'/'+vcf_file, 'r')
	print("{}\{}".format(count, len(vcf_files)))
	for line in vcf.readlines()[1:]:
		variant = break_down_mutation(line)
		#print("{}\t{}".format([i for i in line.split('\t') if i in ['SNP','DEL','INS','LSP']], str(variant)))
	count += 1

def check_variant_commercial(variant):
	drug, gene_name, codon_position, type_change = variant.drug, variant.gene_name, variant.codon_position, variant.type_change
	return_variable = False
	if(drug == 'INH' and gene_name == 'katG' and codon_position == 315):
		return_variable = True
	elif(drug == 'INH' and 'inhA' in gene_name and codon_position in [-15, -16,-8]):
		return_variable = True
	elif(drug == 'RIF'  and gene_name == 'rpoB' and codon_position in list(range(505-81,534-81))):
		return_variable = True
	elif(drug in ['KAN', 'AMK','CAP', 'SLIS'] and gene_name == 'rrs' and codon_position in [1401,1402]):
		return_variable = True
	elif(drug in ['KAN','AMK','CAP', 'SLIS'] and gene_name == 'eis' and codon_position in [-10,-11,-12,-13,-14,-37]):
		return_variable = True
	#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
	elif(drug in ['LEVO', 'CIP','OFLX', 'FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
		return_variable = True
	#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
	elif(drug in ['LEVO', 'CIP','OFLX', 'FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
		return_variable = True

	if(return_variable and type_change[0] == type_change[1]):
		return True, 'synonymous'
	elif(return_variable and type_change[0] != type_change[1]):
		return False, 'asynonymous'
	else:
		return False

# def produce_variant(line, search = 0, search_drug = None):

# 	if(search):
# 		line = line.split('\t')
# 		ID = line[5]
		
# 		#We do some sneaky stuff here -- if & we just split it up and recursively call this function on each part
# 		#In effect we create two lines from one
# 		if('&' in ID):
# 			#print("OOOOHHH you got the speeccciialll caaassseeee")
# 			multiple = ID.split('&')
# 			for i in multiple:
# 				copy = line[:]
# 				copy[5] = i
# 				if(produce_variant('\t'.join(copy), 1, search_drug)):
# 					return True	
# 			return False
	
# 		type_change_info = ID.split('_') 
		
# 		#print(line)
# 		#print(type_change_info)
#         	if(type_change_info[0] == 'SNP' and type_change_info[1] == 'CN' and 'LSP' not in type_change_info):
#                 	gene_name, codonAA = type_change_info[5], type_change_info[4]
# 			if('-' in codonAA):
# 				codon_position = [i for i in re.findall('\d*', codonAA) if i]
# 				type_change = ''
# 			else:
#                 		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
#         	elif(type_change_info[0] == 'SNP' and type_change_info[1] == 'P'):
#                 	gene_name, codonAA = type_change_info[4], type_change_info[3]
#                 	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
#         	elif(type_change_info[0] == 'SNP' and type_change_info[1] == 'N'):
#                 	gene_name, codonAA = type_change_info[4], type_change_info[3]
#                 	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
#         	elif(type_change_info[0] == 'SNP' and type_change_info[1] == 'I' and 'DEL' not in type_change_info and 'INS' not in type_change_info):
#                 	gene_name, codonAA = type_change_info[4], type_change_info[3]
# 			if('-' in codonAA):
# 				codon_position = [i for i in re.findall('\d*', codonAA) if i]
# 				type_change = ''
# 			else:
#                 		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
#         	elif(type_change_info[0] == 'INS' and type_change_info == 'F'):
#                 	gene_name, codon_position = type_change_info[4], type_change_info[3]
#                 	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','')
#         	elif(type_change_info[0] == 'DEL' and type_change_info[1] == 'F'):
#                 	gene_name, codon_position = type_change_info[4], type_change_info[3]
#                 	type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','')
#         	else:
# 				#So this is a non-standard type, but it may still be relevant so we can just use regex
# 				#for gene name we can always assume its in the same position
# 				gene_name = type_change_info[4]
# 				#However, for the codon position we use regex any numbers will be codon position
# 				codon_position = [i for i in re.findall('\d*', type_change_info[3]) if i][0]
# 				type_change = ''
# 	        		#print("THERES SOMETHING ELSE LURKING HERE....RUN! {}".format(line))
			
# 			#return False
# 		#print(line)
# 		#print(type_change_info)
# 		#if(len(codon_position) > 30 or '&' in ID):
# 			#print("CODON MESSED UP: {}\t{}".format(line,codon_position))
# 			#return False

# 	if(not search):
# 		line = line.replace('\n','')
# 		data = line.split('_')
# 		drug, type_change_info = data[0], data[1:]

# 		if('_u' in line):
# 			if('promoter' in line and '.' in line):
# 				# print(line)
# 				type_change_info = ['SNP'] + ['P'] + ['_'] + ['_'] + type_change_info
# 				type_change_info[3] = type_change_info[5]
# 				# print(type_change_info)
# 			else:
# 				# print(line)
# 				type_change_info = ['SNP'] + ['CN'] + ['_'] + ['_'] + type_change_info
# 				# print(type_change_info)
# 		# type_change_info = type_resistance.split('_')
# 		if(type_change_info[1] == 'CN'):
# 			gene_name, codonAA = type_change_info[4], type_change_info[5]
# 			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
# 		elif(type_change_info[1] == 'P'):
# 			gene_name, codon_position = type_change_info[4], type_change_info[3]
# 			if('.' in codon_position):
# 				split = codon_position.split('.')
# 				codon_position = split[1]
# 				type_change = split[0]
# 			if('.' in gene_name):
# 				gene_name = 'promoter-'+gene_name.replace('.','-')
# 			else:
# 				gene_name = 'promoter-'+gene_name
# 		elif(type_change_info[1] == 'N'):
# 			gene_name, codonAA = type_change_info[4], type_change_info[3]
# 			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
# 			#type_change = None
# 		elif(type_change_info[1] == 'I'):
# 			gene_name, codon_position = type_change_info[4], type_change_info[3]
# 			if('.' in codon_position):
# 	                        split = codon_position.split('.')
# 	                        codon_position = split[1]
# 	                        type_change = split[0]
# 			if('.' in gene_name):
# 				gene_name = 'inter-'+gene_name.replace('.','-')
# 		elif(type_change_info[0] == 'INS' and type_change_info[1] == 'F'):
# 			gene_name, codonAA = type_change_info[4], type_change_info[3]
# 			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
# 			#type_change = None
# 		elif(type_change_info[0] == 'DEL' and type_change_info[1] == 'F'):
# 			gene_name, codonAA = type_change_info[4], type_change_info[3]
# 			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
# 			#type_change = None
# 		else:
# 			print("THERES SOMETHING ELSE LURKING HERE....RUN! {}".format(line))

# 	"""

# 	From presentation
# 	INH katG codon 315
# 	INH inhA promoter -15, -16, -8
# 	RIF rpoB codon 505-534
# 	AMINOGLYCOSIDES rrs 1401, 1402
# 	AMINOGLYCOSIDES eis promoter -10 to -14, -37
# 	FLQ gyrA 88-94
# 	FLQ gyrB 500-541
# 	"""
# 	#THIS IS WHERE the info for commercial is stored -- we check if commercial!
# 	def check(drug, gene_name, codon_position, search, type_change):
# 		add = False
# 		if(drug == 'INH' and gene_name == 'katG' and codon_position == 315):
# 			drug = 'INH'
# 			add = True
# 		elif(drug == 'INH' and 'inhA' in gene_name and codon_position in [-15, -16,-8]):
# 			drug = 'INH'
# 			add = True
# 		elif(drug == 'RIF'  and gene_name == 'rpoB' and codon_position in list(range(505-81,534-81))):
# 			drug = 'RIF'
# 			add = True
# 		elif(drug in ['KAN', 'AMK','CAP', 'SLIS'] and gene_name == 'rrs' and codon_position in [1401,1402]):
# 			drug = 'SLIS'
# 			add = True
# 		elif(drug in ['KAN','AMK','CAP', 'SLIS'] and gene_name == 'eis' and codon_position in [-10,-11,-12,-13,-14,-37]):
# 			drug = 'SLIS'
# 			add = True
# 		#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
# 		elif(drug in ['LEVO', 'CIP','OFLX', 'FLQ'] and gene_name == 'gyrA' and codon_position in [88,89,90,91,92,93,94]):
# 			drug = 'FLQ'
# 			add = True
# 		#elif(drug in ['LEVO','FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
# 		elif(drug in ['LEVO', 'CIP','OFLX', 'FLQ'] and gene_name == 'gyrB' and codon_position in list(range(500,541))):
# 			drug = 'FLQ'
# 			add = True

# 		if(add and search):
# 			print("Checking: {}".format(type_change))
# 			if(type_change[0] == type_change[1]):
# 				print("BOIIIII DAMN WE FOUND IT A SYNONYMOUS AND WE WERE ABOUT TO APPROVVEEE")
# 				print(type_change)
# 			return True
# 		elif(not add and search):
# 			return False
# 		return add


# 	if(search and search_drug):
# 		drug = search_drug

# 	add = False
# 	codon_positionn = None
# 	if(isinstance(codon_position, list)):
# 		codon_positionn = [check(drug, gene_name, int(i), search, type_change) for i in codon_position]
# 		if False in codon_positionn:
# 			return False
# 		else:
# 			return True
# 	else:
# 		codon_position = int(codon_position)
# 		add = check(drug, gene_name, codon_position, search, type_change)
# 		if(not add):
# 			return False
# 		elif(add and search):
# 			return True

# 	if(add):
# 		if(drug in ['LEVO','CIP','OFLX']):
# 		#if(drug in ['LEVO']):
# 			drug = 'FLQ'
# 		elif(drug in ['KAN','AMK','CAP']):
# 			drug = 'SLIS'
# 		if(drug not in drug_to_variants):
# 			drug_to_variants[drug] = []
# 		else:
# 			drug_to_variants[drug].append(line)


# #Complex logic here -- basically we only add to drug_to_variants if its a commercial one that we already found
# for line in [i for i in list(res.columns)[1:] if 'unknown' not in i and 'embA' not in i and 'PAS' not in i and 'PZA' not in i]:
# 	produce_variant(line)


# variants_unique_WGS = set()
# strains_analyzed = set()
# commercial
# WGS = 


# for d in ['INH','RIF','SLIS','FLQ']:
# #for d in ['FLQ']:
# 	"""
# 		FLQ includes CIP, OFLX, LEVO
# 		SLIS includes KAN, AMK, CAP

# 	"""
# 	if(d not in drug_to_variants):
# 		print("{} not detected at all".format(d))
# 	else:
# 		mutations_interest = drug_to_variants[d]
# 		print(mutations_interest)
# 		if(d == 'INH'):
# 			d = 'ISONIAZID'
# 		sub = combined[combined[d] == 1].copy()
# 		#sub = combined[combined['LEVOFLOXACIN'] == 1].copy()
# 		total = 0

# 		total_count   = len(sub.index)
# 		count = 0
# 		for index, row in sub.iterrows():
# 			det = False
# 			strains_analyzed.add(row['strain'])

# 			def perform_commercial_test(row, mutations_interest):
# 				#Step 1: Search for commercial we already found (i.e. already in list, we already looked for it)
# 				for mutation in mutations_interest:
# 					if(int(row[str(mutation)]) == 1):
# 						det = True
# 						file.write('{}\t{}\t{}\n'.format(row['strain'], str(mutation), 'commercial'))
# 						break

# 				#Step 2: If we found no commercial, then search the vcf for a mutation in the genes we didn't already find but are tested by commercial drugs
# 				if(not det):
# 					#Go through lilne by line in vcf
# 					vcf_file = [i for i in vcf_files if row['strain']+'.var'== i]
# 					if(len(vcf_file) > 1):
# 						raise Exception('More than one vcf file found for {}'.format(row['strain']))
# 					elif(not vcf_file):
# 						raise Exception('No vcf file found for {}'.format(row['strain']))
					
# 					name = vcf_file[0]
# 					vcf = open(vcf_directory+'/'+name, 'r')
# 					for line in vcf.readlines()[1:]:
# 						if(produce_variant(line, 1, d)):
# 							print(line)
# 							print("wow good thing we did this found a commercial!")
# 							det = True
# 							file.write('{}\t{}\t{}\n'.format(row['strain'], line.split('\t')[5], 'commercial'))
# 							break

# 			def perform_WGS_test(row):
# 				if(d == 'ISONIAZID'):
# 					d = 'INH'
# 				WGS_only = []
# 				for i in list(sub.columns):
# 					if i not in mutations_interest:
# 						if 'unknown' not in i:
# 							if '_' in i:
# 								if(d == 'INH' or d == 'RIF'):
# 									if(d in i):
# 										WGS_only.append(i)
# 								else:
# 									if(d == 'FLQ'):
# 										if('CIP' in i or 'OFLX' in i or 'LEVO' in i):
# 										#if('LEVO' in i):
# 											WGS_only.append(i)
# 									elif(d == 'SLIS'): 
# 										if('KAN' in i or 'AMK' in i or 'CAP' in i):
# 											WGS_only.append(i)

# 				for mutation in WGS_only:
# 					if(int(row[str(mutation)]) == 1):
# 						variants_unique_WGS.add(row['strain'])
# 						det = True
# 						file.write('{}\t{}\t{}\n'.format(row['strain'], str(mutation), 'WGS'))
# 						break

# 			#Step 4: If after all that commercial and WGS can't find it then the strain is undetected with current information
# 			if(not det):
# 				undetected += 1

# 			count += 1
# 			print("{}\t{}".format(count, total_count))
# 			total += 1
# 		total = float(total)
# 		print("DRUG: {} \t COMMERCIAL DETECTED: {} \t TOTAL: {} \t PERCENT DETECTED: {}".format(d, detected_commerical, total, detected_commerical/total))
# 		print("DRUG: {} \t WGS DETECTED: {} \t TOTAL: {} \t PERCENT DETECTED: {}".format(d, detected_WGS, total, detected_WGS/total))
# 		print("DRUG: {} \t UNDETECTED: {} \t TOTAL: {} \t PERCENT UNDETECTED: {}".format(d, undetected, total, undetected/total))


# print("TOTAL NUMBER ISOLATES RESISTANCE IDENTIFIED VIA WGS: {} \ {}".format(len(list(variants_unique_WGS)),len(list(strains_analyzed))))





