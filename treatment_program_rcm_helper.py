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
			return_string = self.AA_change+self.codon_location
		else:
			return_string = self.codon_location
		return "{}_{}".format(self.gene_name, return_string)

def annotate(df):
	"""
	These next definitions annotate the strain_info.tsv file with drug categories we use later on in analysis

	Conditions for various categories 
	1 ==> We have the data that shows that this isolate should be in this category
	0 ==> We have the data that shows that this isolate should NOT be in this category
	-1 ==> We don't have enough data to say anything :(

	"""

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

	df['FQ'] = df.apply(lambda row: FLUOROQUINOLONES(row),axis=1)

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

	#To be MDR the isolate has to be resistant to both Isoniazid and RIF
	def MDR(row):
		if(row['RIF'] == 1 and row['ISONIAZID'] == 1):
			return 1
		elif(row['RIF'] == -1 or row['ISONIAZID'] == -1):
			return -1
		else:
			return 0

	df['MDR'] = df.apply(lambda row: MDR(row),axis=1)

	#To be XDR the isolate has to be MDR and resistant to Fluoroquinolones and SLIS
	def XDR(row):
		if(row['MDR'] == 1 and row['FQ'] == 1 and row['SLIS'] == 1):
			return 1
		elif(row['MDR'] == -1 or row['FQ'] == -1 or row['SLIS'] == -1):
			return -1
		else:
			return 0

	df['XDR'] = df.apply(lambda row: XDR(row),axis=1)

	#For these next categories we only look at isolates that have data on isoniazid and rifamycin (i.e df[(df['ISONIAZID'] > -1 ) & (df['RIFAMYCIN'] > 1)])
	#To be susceptible the isolate has to have data on isoniazid and rif and then not have any data that shows it is resistant
	def S(row):
		if(row['ISONIAZID'] > -1 and row['RIF'] > -1):
			if(row['ISONIAZID'] < 1 and row['RIF'] < 1 and row['FQ'] < 1 and row['SLIS'] < 1 and row['PAS'] < 1 and row['CYS'] < 1 and row['ETH'] < 1 and row['ETHAMBUTOL'] < 1 and row['STREPTOMYCIN'] < 1 and row['PYRAZINAMIDE'] < 1):
				return 1
			else:
				return 0
		else:
			return -1

	df['S'] = df.apply(lambda row: S(row),axis=1)

	#For Isoniazid mono resistant we look at isolates which have data on isoniazd and rif -- if the isolate has this data we then make sure it is only resistant to INH and susceptible to rif 
	#and susceptible to any other drug the isolate was tested for except strep where we allow it to be resistant to both INH and strep and still be INH mono
	def INH_MONO(row):
		if(row['ISONIAZID'] > -1 and row['RIF'] > -1):
			if(row['ISONIAZID'] == 1 and row['RIF'] == 0 and row['FQ'] < 1 and row['SLIS'] < 1 and row['PAS'] < 1 and row['CYS'] < 1 and row['ETH'] < 1 and row['ETHAMBUTOL'] < 1 and row['STREPTOMYCIN'] <= 1 and row['PYRAZINAMIDE'] < 1):
				return 1
			else:
				return 0
		else:
			return -1

	df['INH_MONO'] = df.apply(lambda row: INH_MONO(row),axis=1)

	#Similar logic to the INH Mono except now we require that is only be resistant to Strep and susceptible to all other drugs
	def STREP_MONO(row):
		if(row['ISONIAZID'] > -1 and row['RIF'] > -1):
			if(row['ISONIAZID'] == 0 and row['RIF'] == 0 and row['FQ'] < 1 and row['SLIS'] < 1 and row['PAS'] < 1 and row['CYS'] < 1 and row['ETH'] < 1 and row['ETHAMBUTOL'] < 1 and row['STREPTOMYCIN'] == 1 and row['PYRAZINAMIDE'] < 1):
				return 1
			else:
				return 0
		else:
			return -1

	df['STR_MONO'] = df.apply(lambda row: STREP_MONO(row),axis=1)

	#Similar logic to the previous 2 except now we check if the isolate is not INH mono, strep mono, xdr, or mdr and it is still resistant to some drug then it is other resistant
	def OTHER_R(row):
		if(row['ISONIAZID'] > -1 and row['RIF'] > -1):
			if((row['INH_MONO'] < 1 and row['STR_MONO'] < 1 and row['XDR'] < 1 and row['MDR'] < 1 ) and (row['ISONIAZID'] == 1 or row['RIF'] == 1 or row['FQ'] == 1 or row['SLIS'] == 1 or row['PAS'] == 1 or row['CYS'] == 1 or row['ETH'] == 1 or row['ETHAMBUTOL'] == 1 or row['STREPTOMYCIN'] == 1 or row['PYRAZINAMIDE'] == 1)):
				return 1
			else:
				return 0
		else:
			return -1

	df['OTHER_R'] = df.apply(lambda row: OTHER_R(row),axis=1)

	#Find if an entry has no resistance data
	def NO_DATA(row):
		if(row['S'] < 0 and row['MDR'] < 0 and row['INH_MONO'] < 0 and row['STR_MONO'] < 0 and row['OTHER_R'] < 0 and row['XDR'] < 0):
			return 1
		return 0

	df['NO_DATA'] = df.apply(lambda row: NO_DATA(row),axis=1)
	return df.copy()


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



