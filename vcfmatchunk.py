import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description='Takes in a list of IDS breaks down the codon locations, gene names, and type of mutation and looks for that in the other inputted vcf file')

parser.add_argument('vcf_files', help='directory with vcf files')
parser.add_argument('ref_file', help='the vcf file with reference IDs to look for')
parser.add_argument('interest_file', help='the vcf file with reference IDs of interest of those to look for')
parser.add_argument('strain_info_file', help='the strain_info_file with resistance data for every strain')
parser.add_argument('lineage_snps_file', help='file with snps that are not rcms but lineage specific to discard')

args = parser.parse_args()
ref_file, vcf_directory, interest_file = open(args.ref_file, 'r'), args.vcf_files, open(args.interest_file, 'r')
snps_of_interest = [i.replace('\n','') for i in interest_file.readlines()]
strain_info = pd.read_csv(args.strain_info_file,sep='\t')

#Added filtering to only use strains we had data on at time of analysis
strains_we_have = [i.rstrip() for i in open('final_all_strains','r').readlines()]
strain_info  = strain_info[strain_info['strain'].isin(strains_we_have)]

strain_info = strain_info.set_index('strain')

drugs = set()

drug_mapping = {
		'AMK':'AMIKACIN',
		'CAP':'CAPREOMYCIN',
		'CIP':'CIPROFLOXACIN',
		'EMB':'ETHAMBUTOL',
		'ETH':'ETHIONAMIDE',
		'INH':'ISONIAZID',
		'KAN':'KANAMYCIN',
		'LEVO':'LEVOFLOXACIN',
		'OFLX':'OFLOXACIN',
		'PAS':'PARA-AMINOSALICYLIC_ACID', 
		'PZA':'PYRAZINAMIDE',
		'RIF':'RIFAMPICIN',
		'STR':'STREPTOMYCIN'
	}

class Variant:

	def __init__(self, gene_name, codon_location, AA_change, name = None, test_name = None, drug = None):
		self.gene_name = gene_name
		self.codon_location = codon_location
		self.AA_change = AA_change
		self.name = name
		if(test_name):
			self.is_top = (drug+'_'+test_name) in snps_of_interest
		else:
			self.is_top = False
		self.drug = drug
		# print(test_name)
	def compare_variant(self, variant):
		return (self.gene_name == variant.gene_name) and (self.codon_location == variant.codon_location) and (self.AA_change == variant.AA_change)

	def __str__(self):
		if(self.AA_change):
			return_string = self.AA_change[0]+self.codon_location+self.AA_change[1]
		else:
			return_string = self.codon_location
		return "{}_{}".format(self.gene_name, return_string)


"""
Example piece of data:

INH	SNP_CN_2155168_CG_katG_S315T
INH	SNP_P_1673425_CT.15_fabG1
STR	SNP_N_1472359_A514C_rrs
PZA	INS_F_2288725_i516C_pncA
PZA	DEL_F_2289069_d172A_pncA_F58L
STR	SNP_I_1473637_A.21_rrs.rrl
"""
drug_to_variants = {}
all_variants_to_compare = []
drugs_to_gene = {}

for line in ref_file.readlines():
	line = line.replace('\n','')
	data = line.split('\t')
	drug, type_resistance = data[0], data[1]
	if(drug not in drug_to_variants):
		drug_to_variants[drug] = []
	type_change_info = type_resistance.split('_')
	if(type_change_info[1] == 'CN'):
		gene_name, codonAA = type_change_info[4], type_change_info[5]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1]
	elif(type_change_info[1] == 'P'):
		gene_name, codon_position = type_change_info[4], type_change_info[3]
		if('.' in codon_position):
			split = codon_position.split('.')
			codon_position = split[1]
			type_change = split[0]
		if('.' in gene_name):
			gene_name = 'promoter-'+gene_name.replace('.','-')
		else:
			gene_name = 'promoter-'+gene_name
	elif(type_change_info[1] == 'N'):
		gene_name, codonAA = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
		#type_change = None
	elif(type_change_info[1] == 'I'):
		gene_name, codon_position = type_change_info[4], type_change_info[3]
		if('.' in codon_position):
                        split = codon_position.split('.')
                        codon_position = split[1]
                        type_change = split[0]
		if('.' in gene_name):
			gene_name = 'inter-'+gene_name.replace('.','-')
	elif(type_change_info[0] == 'INS' and type_change_info[1] == 'F'):
		gene_name, codon_position = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
		#type_change = None
	elif(type_change_info[0] == 'DEL' and type_change_info[1] == 'F'):
		gene_name, codon_position = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
		#type_change = None
	else:
		print("THERES SOMETHING ELSE LURKING HERE....RUN! {}".format(line))

	drugs.add(drug)
	variant  = Variant(gene_name, codon_position, type_change, drug+'_'+type_resistance.replace('\n',''), type_resistance.replace('\n',''), drug)
	if(drug in drugs_to_gene):
		drugs_to_gene[drug].add(gene_name)
	else:
		drugs_to_gene[drug] = set([gene_name])

	drug_to_variants[drug].append(variant)
	all_variants_to_compare.append(variant)

###ADD LEVO gyrB gene###
drugs_to_gene['LEVO'].add('gyrB')

####ADD PAS inter-thyX-hsdS.1 and folC GENE###
drugs_to_gene['PAS'].add('folC')
drugs_to_gene['PAS'].add('inter-thyX-hsdS.1')

"""
Now break down lineage specific snps

"""

lineage_snps = []

for line in open(args.lineage_snps_file,'r').readlines():
	line = line.replace('\n','')
	type_change_info = line.split('_')
	if(type_change_info[1] == 'CN'):
		gene_name, codonAA = type_change_info[5], type_change_info[4]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
	elif(type_change_info[1] == 'P'):
		gene_name, codonAA = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
	elif(type_change_info[1] == 'N'):
		gene_name, codonAA = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
	elif(type_change_info[1] == 'I'):
		gene_name, codonAA = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
	elif(type_change_info[0] == 'INS' and type_change_info == 'F'):
		gene_name, codon_position = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
	elif(type_change_info[0] == 'DEL' and type_change_info[1] == 'F'):
		gene_name, codon_position = type_change_info[4], type_change_info[3]
		type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 

	test = Variant(gene_name, codon_position, type_change)
	lineage_snps.append(test)


""" SNP_CN_7585_G284C_S95T_gyrA

SNP_P_781395_T165C_promoter-rpsL

SNP_N_1474571_G914A_rrl

INS_CF_2289050_i192A_65I_pncA

"""

drugs = list(drugs)
print(drugs_to_gene)
master = {}
for name in [i for i in os.listdir(vcf_directory) if '.var' in i]:
	vcf_mapping = {}
	vcf_file = open(vcf_directory+'/'+name, 'r')
	all_variants_identified = []

	for line in vcf_file.readlines()[1:]:
		line = line.split('\t')
		ID = line[5]
		type_change_info = ID.split('_')
		if(type_change_info[1] == 'CN'):
			gene_name, codonAA = type_change_info[5], type_change_info[4]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
		elif(type_change_info[1] == 'P'):
			gene_name, codonAA = type_change_info[4], type_change_info[3]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
		elif(type_change_info[1] == 'N'):
			gene_name, codonAA = type_change_info[4], type_change_info[3]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
		elif(type_change_info[1] == 'I'):
			gene_name, codonAA = type_change_info[4], type_change_info[3]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1] 
		elif(type_change_info[0] == 'INS' and type_change_info == 'F'):
			gene_name, codon_position = type_change_info[4], type_change_info[3]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
		elif(type_change_info[0] == 'DEL' and type_change_info[1] == 'F'):
			gene_name, codon_position = type_change_info[4], type_change_info[3]
			type_change,codon_position = codonAA[0]+codonAA[len(codonAA)-1], codonAA[1:len(codonAA)-1].replace('\n','') 
		if('Rv' not in gene_name):
			#None of the gene names has Rv in them so if this shows up...bye!
			#if(gene_name == 'rrs' and codon_position == '1401' and type_change == 'AG'):
			#	print("HOLY ITS RIGHT HERE TOO SO ITS IN DEEEEEEEP")
			test = Variant(gene_name, codon_position, type_change)
			all_variants_identified.append(test)

	drugs_found = set()

	for compare_variant in all_variants_to_compare:
		for variant in all_variants_identified:
			#if(variant.gene_name == 'rrs' and variant.codon_location == '1401' and variant.AA_change == 'AG'):
			#	print("HOLY THIS IS IT {}\t{}\t{}".format(compare_variant.gene_name, compare_variant.codon_location, compare_variant.AA_change))	 
			if(compare_variant.compare_variant(variant)):
				if(compare_variant.is_top):
					vcf_mapping[compare_variant.name] = 1
					drugs_found.add(compare_variant.drug)
				else:
					vcf_mapping['{}_other_snp'.format(compare_variant.drug)] =  1
					drugs_found.add(compare_variant.drug)
				print("{} we found this {}".format(name, compare_variant.name))
				break
			else:
				if(compare_variant.is_top):
					vcf_mapping[compare_variant.name] = 0
				else:
					if('{}_other_snp'.format(compare_variant.drug) not in vcf_mapping):
						vcf_mapping['{}_other_snp'.format(compare_variant.drug)] = 0
	name = name.replace('.var','')

	#Now we go after the unknowns -- checking which drugs have unknown and if they are unknown 
	for drug in drugs:
		#Check if strain exists in table of resistant to strain and if it does is it resistant?
		check_if_exists = name  in strain_info.index
		is_resistant = False
		if(check_if_exists):
			real_drug = drug_mapping[drug]
			if(int(strain_info.at[name,real_drug]) > 0):
				is_resistant = True

		#Okay now if we haven't found a RCM for that drug 
		if(drug not in drugs_found):
			# print("secretstuff")
			#If the isolate is resistant for that drug then we have an unknown
			if(is_resistant):
				#Now cycle through all the variant we found for this isolate and if they match the drug and associate gene location we add it as another snp
				#print("OOOHHHH WE CHECKIN")
				explained = False
				for variant in all_variants_identified:
					if( variant.gene_name in drugs_to_gene[drug] ):
						lineage_mutation = False
						for lineage_variant in lineage_snps:
							if(lineage_variant.compare_variant(variant)):
								lineage_mutation = True
								print("OOOPS so close but lineage mutation")
						if(not lineage_mutation):
							#OOOYEEE we found one -- create a new category for this
							print("WOH WE FOUND A NEW SNPPPPPPPP")
							vcf_mapping[drug+'_'+str(variant)+'_u'] = 1
							explained = True
				if(not explained):
					#So if even after searching all the variants we were not able to find a new snp then this is truly unknown
					vcf_mapping['{}_unknown'.format(drug)] = 1	
			else:
				#If the isolate is not resistant for that drug and we found nothing, we simply don't have enough data to say anything (i.e. there can be a row of all 0s)
				vcf_mapping['{}_unknown'.format(drug)] = 0
		else:
			#If we did find a RCM for that drug then this isolate cannot be categorized as unknown for that drug
			vcf_mapping['{}_unknown'.format(drug)] = 0

	master[name] = vcf_mapping


results = pd.DataFrame.from_dict(master).T.copy()

#change all NaNs to zero cause we didn't find those in those isolates
results = results.fillna(0)

results.to_csv('results_modified_unknown',sep='\t')
