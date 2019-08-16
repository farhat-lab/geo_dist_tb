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

	return Variant(gene_name, codon_position, type_change, drug=drug)  

def check_variant_commercial(variant, drug):
	gene_name, codon_position, type_change = variant.gene_name, variant.codon_location, variant.AA_change
	if(drug == 'ISONIAZID'):
		drug = 'INH'
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
		return 'synonymous'
	elif(return_variable and type_change[0] != type_change[1]):
		return 'asynonymous'
	else:
		return False

file_resultat = open('commercial_results','w')
file_resultat.write('drug\tmutation\tresistant\tsusceptible\tsynonymous\tasynonymous\n')
for drug in ['ISONIAZID','RIF','SLIS','FQ']:
	#memoization
	memoization = {}
	count = 0
	#Only look at strains with data
	sub = combined[(combined['NO_DATA'] != 1) &(combined[drug] != -1)]
	for index, row in sub.iterrows():
		print("{}\t{}".format(count, len(list(sub.iterrows()))))
		vcf = open(vcf_directory+'/'+row['strain']+'.var','r')
		for line in vcf.readlines()[1:]:
			mutation = [i for i in line.split('\t') if 'SNP' in i or 'INS' in i or 'DEL' in i or 'LSP' in i][0]
			if(mutation in memoization):
				if(memoization[mutation]):
					file_resultat.write(memoization[mutation])
			else:
				commercial = check_variant_commercial(break_down_mutation(line), drug)
				if(commercial):
					if(row[drug] != 1):
						if(commercial == 'synonymous'):
							memoization[mutation] = "{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 0,1,1,0)
							print("{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 0,1,1,0))
							file_resultat.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 0,1,1,0))
						elif(commercial == 'asynonymous'):
							memoization[mutation] = "{}\t{}\t{}\t{}\t{}\t{\n}".format(drug, line.split('\t')[5], 0,1,0,1)
							print("{}\t{}\t{}\t{}\t{}\t{\n}".format(drug, line.split('\t')[5], 0,1,0,1))
							file_resultat.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 0,1,0,1))
					else:
						if(commercial == 'synonymous'):
							memoization[mutation] = "{}\t{}\t{}\t{}\t{}\t{\n}".format(drug, line.split('\t')[5], 1,0,1,0)
							print("{}\t{}\t{}\t{}\t{}\t{\n}".format(drug, line.split('\t')[5], 1,0,1,0))
							file_resultat.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 1,0,1,0))
						elif(commercial == 'asynonymous'):
							memoization[mutation] = "{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 1,0,0,1)
							print("{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 1,0,0,1))
							file_resultat.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(drug, line.split('\t')[5], 1,0,0,1))
				else:
					memoization[mutation] = False
		count += 1






