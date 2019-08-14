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

def annotate(df):
	#Rifamycin category requires that isolate be resistant to either rifabutin or rifampicin
	def RIFAMYCIN(row):
		if(row['RIFABUTIN'] == 1 or row['RIFAMPICIN'] == 1):
			return 1
		elif(row['RIFABUTIN'] == 0 or row['RIFAMPICIN'] == 0):
			return 0
		else:
			return -1

	df['RIF'] = df.apply(lambda row: RIFAMYCIN(row),axis=1)

	#Fluoroquinolones category required that isolate be resistant to either moxifloxacin, ciprofloxacin, or ofloxacin
	def FLUOROQUINOLONES(row):
		if(row['MOXIFLOXACIN'] == 1 or row['CIPROFLOXACIN'] == 1 or row['OFLOXACIN'] == 1):
			return 1
		elif(row['MOXIFLOXACIN'] == 0 or row['CIPROFLOXACIN'] == 0 or row['OFLOXACIN'] == 0):
			return 0
		else:
			return -1

	df['FLQ'] = df.apply(lambda row: FLUOROQUINOLONES(row),axis=1)

	#Slis category requires that isolate be resistant to either kanamycin or amikacin or capreomycin
	def SLIS(row):
		if(row['KANAMYCIN'] == 1 or row['AMIKACIN'] == 1 or row['CAPREOMYCIN'] == 1):
			return 1
		elif(row['KANAMYCIN'] == 0 or row['AMIKACIN'] == 0 or row['CAPREOMYCIN'] == 0):
			return 0
		else:
			return -1

	df['SLIS'] = df.apply(lambda row: SLIS(row),axis=1)

	#Ethionamide category requires that isolates be resistant to Ethionamide or prothionamide
	def ETHIONAMIDE(row):
		if( (row['ETHIONAMIDE'] == 1) or (row['PROTHIONAMIDE'] == 1) ):
			return 1
		elif( (row['ETHIONAMIDE'] == 0) or (row['PROTHIONAMIDE'] == 0) ):
			return 0
		else:
			return -1

	df['ETH'] = df.apply(lambda row: ETHIONAMIDE(row),axis=1)

	#Cycloserine category requires that isolates be resistant to cycloserine
	def CYCLOSERINE(row):
		if( (row['CYCLOSERINE'] == 1) ):
			return 1
		elif( (row['CYCLOSERINE'] == 0) ):
			return 0
		else:
			return -1

	df['CYS'] = df.apply(lambda row: CYCLOSERINE(row),axis=1)

	#Para-aminosalicyclic acid category requires that isolates be resistant to para-aminosalicylic acid
	def PARA_AMINOSALICYLIC_ACID(row):
		if( row['PARA-AMINOSALICYLIC_ACID'] == 1):
			return 1
		elif( row['PARA-AMINOSALICYLIC_ACID'] == 0 ):
			return 0
		else:
			return -1

	df['PAS'] = df.apply(lambda row: PARA_AMINOSALICYLIC_ACID(row),axis=1)



drugs = [ 'AMK', 'CAP', 'CIP', 'EMB', 'INH', 'KAN', 'LEVO', 'OFLX', 'PAS', 'PZA', 'RIF', 'STR']

drug_mapping = {
		'AMK':'AMIKACIN',
		'CAP':'CAPREOMYCIN',
		'CIP':'CIPROFLOXACIN',
		'EMB':'ETHAMBUTOL',
		'INH':'ISONIAZID',
		'KAN':'KANAMYCIN',
		'LEVO':'LEVOFLOXACIN',
		'OFLX':'OFLOXACIN',
		'PAS':'PARA-AMINOSALICYLIC_ACID', 
		'PZA':'PYRAZINAMIDE',
		'RIF':'RIFAMPICIN',
		'STR':'STREPTOMYCIN',
		'ETH':'ETHIONAMIDE'
	}



