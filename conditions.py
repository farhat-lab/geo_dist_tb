import pandas as pd 
import numpy as np  
import statsmodels.stats.proportion as sp
import operator
import scipy.stats as sc
import itertools
import statsmodels.api as sm
import matplotlib.pyplot as plt

df = pd.read_csv('strain_info.tsv',sep='\t')
strains_we_have = [i.rstrip() for i in open('final_all_strains','r').readlines()]
df = df[df['strain'].isin(strains_we_have)]
resistance_mutation = pd.read_csv('results_modified_unknown',sep='\t')
resistance_mutation = resistance_mutation[resistance_mutation['strain'].isin(strains_we_have)].copy()
resistance_mutation = pd.merge(df, resistance_mutation, on='strain', how='inner')
name_with_date = pd.read_csv('strain_name_with_date_unprocessed',sep='\t')[['strain','drug','date']].copy()
name_with_date['drug']=name_with_date['drug'].replace({'MOXIFLOXACIN':'FLQ','CIPROFLOXACIN':'FLQ','OFLOXACIN':'FLQ'})
name_with_date['drug']=name_with_date['drug'].replace({'KANAMYCIN':'SLIs','AMIKACIN':'SLIs','CAPREOMYCIN':'SLIs'})
name_with_date['drug']=name_with_date['drug'].replace({'RIFAMPICIN':'RIF','RIFABUTIN':'RIF','RIFAMPICIN/RIFABUTIN':'RIF'})
name_with_date['drug']=name_with_date['drug'].replace({'ETHIONAMIDE':'ETH','PROTHIONAMIDE':'ETH'})
name_with_date['drug']=name_with_date['drug'].replace({'CYCLOSERINE':'CYS'})
name_with_date['drug']=name_with_date['drug'].replace({'ETHAMBUTOL':'EMB'})
name_with_date['drug']=name_with_date['drug'].replace({'PARA-AMINOSALICYLIC_ACID':'PAS'})
name_with_date['drug']=name_with_date['drug'].replace({'PYRAZINAMIDE':'PZA'})
name_with_date['drug']=name_with_date['drug'].replace({'STREPTOMYCIN':'STR'})
name_with_date['drug']=name_with_date['drug'].replace({'ISONIAZID':'INH'})


##########Support Functions and Definitions#########

def confidence_interval_lower(num, dem):
	bounds = sp.proportion_confint(num, dem, method = 'beta')
	return bounds[0]

def confidence_interval_upper(num, dem):
	bounds = sp.proportion_confint(num, dem, method = 'beta')
	return bounds[1]

def cal(num, dem):
	low = confidence_interval_lower(num, dem)
	up = confidence_interval_upper(num, dem)
	return '{}-{}'.format(low, up)

country_to_region = {
	'Europe':['Georgia','Romania','Belarus','Sweden','Ireland','Denmark','Spain','Netherlands','Russia','Moldova','Portugal',
	'United Kingdom','Switzerland','Germany','Estonia'],
	'Asia':['Thailand','China','Philippines','Nepal','Kazakhstan','Uzbekistan','Pakistan','Burma','India','Indonesia',
	'Iran','South Korea','Turkmenistan','Azerbaijan','Bangladesh','Vietnam'],
	'Africa':['Djibouti','Gambia','Mali','Democratic Republic of the Congo','Morocco','Guinea','Rwanda','Dominican Republic',
	'Ghana','South Africa','Nigeria','United Repubic of Tanzania','Uganda','Sierra Leone','Swaziland',
	'Malawi'],
	'North America':['Canada','United States of America'],
	'South America':['Peru','Colombia','Brazil']
}

WHO_incidence = {
	'Belarus':[0.71,.67,.75],
	'Canada': [0.01,0.01,0.02],
	'China':[0.08,0.07,0.09],
	'Germany':[0.03,0.02,0.05],
	'India':[0.05,0.04,0.06],
	'Iran':[0.02,0.01,0.02],
	'Malawi':[0.01,0.01,0.02],
	'Mali':[0.03,0.02,0.04],
	'Moldova':[0.34,0.31,0.36],
	'Netherlands':[0.02,0.001,0.01],
	'Peru':[0.09,0.09,0.1],
	'Romania':[0.05,0.05,0.06],
	'Russia':[0.43,0.43,0.44],
	'Sierra Leone':[0.03,0.02,0.04],
	'South Africa':[0.04,0.04,0.05],
	'South Korea':[0.04,0.03,0.04],
	'Swaziland':[0.1,0.09,0.11],
	'Thailand':[0.03,0.03,0.04],
	'Turkmenistan':[0.22,0.19,0.24],
	'Uganda':[0.02,0.02,0.03],
	'United Kingdom':[0.02,0.01,0.02],
	'United States of America':[0.02,0.01,0.02],
	'Uzbekistan':[0.33,0.29,0.34]
}

date_drug_introduced = {
	#IN years ago
	'INH':2019-1952,
	'RIF':2019-1965,
	'EMB':2019-1965,
	'STR':2019-1944,
	'PZA':2019-1952,
	'FLQ':2019-1985,
	'SLIs':2019-1958,
	'ETH':2019-1960,
	'CYS':2019-1955
}

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


"""
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

#######################################################################################################

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




###################################

def figure_one():
	number_started_off_with = df['strain'].nunique()
	number_strain_country_info = df[df['country'] != 'Not Provided']['strain'].nunique()
	number_passed_WGS_criteria = resistance_mutation['strain'].nunique()
	number_original_country = df[df['country'] != 'Not Provided']['country'].nunique()
	number_of_resistant_strains_dated = name_with_date['strain'].nunique()
	number_of_dates = name_with_date['strain'].count()


	name_with_date_unclean = pd.read_csv('strain_name_with_date_unprocessed',sep='\t').copy()
	number_of_lineage_drug_groups = (name_with_date_unclean['country']+name_with_date_unclean['lineage']+name_with_date_unclean['drug']).nunique()
	number_of_drugs_groups_dated = name_with_date['drug'].nunique() 

	lineages = df[df['country'] != 'Not Provided'].copy()

	lineages['SUM'] = 1
	test = lineages.replace({-1:0}).groupby('country').sum().reset_index().copy()
	countries_include =  list(set(test[test['SUM'] >= 10]['country'].values))
	lineages = lineages[lineages['country'].isin(countries_include)].copy()



	print("Number we started off with {}".format(number_started_off_with))
	print("Number with country information {}".format(number_strain_country_info))
	print("Number countries in data before any filtering {}".format(number_original_country))
	print("Number that did not fail WGS QC criteria {}".format(number_passed_WGS_criteria))
	print("Number of resistant strains we analyzed for dating {}".format(number_of_resistant_strains_dated))
	print("Number of dates {}".format(number_of_dates))
	print("Number of groups created {}".format(number_of_lineage_drug_groups))
	print("Number of drugs dated within groups (in groups defined in methods) {}".format(number_of_drugs_groups_dated - 1)) ##Note -1 is because we don't include PAS
	print("Number used for lineage analysis in countries more 10 isolates {}".format(lineages['strain'].nunique()))

	lineages = df[df['country'] != 'Not Provided'].copy()

	lineages['SUM'] = 1
	lineages = lineages[lineages['NO_DATA'] != 1]
	test = lineages.replace({-1:0}).groupby('country').sum().reset_index().copy()
	countries_include =  list(set(test[test['SUM'] >= 10]['country'].values))
	lineages = lineages[lineages['country'].isin(countries_include)].copy()


	print("Number used in resistance analysis with INH/RIF data {}".format(lineages['strain'].nunique()))


# figure_one()

##CALCULATE NUMBER COUNTRIES IN SET##
def return_number_countries():
	print(len(list(set(df['country'].values))))

##CALCULATE NUMBER OF ISOLATES##
def return_number_isolates():
	print(df['strain'].nunique())
	print(df[df['NO_LINEAGE'] == 1]['strain'].nunique())
	print(df[df['country'] != 'Not Provided']['strain'].nunique())


###LINEAGE ANALYSIS###
def generate_lineage_data():
	lineages = df.copy()
	lineages['SUM'] = 1
	lineages = lineages.replace({-1:0}).groupby('country').sum()
	lineages[['lineage1','lineage2','lineage3','lineage4','lineage5','lineage6','lineage7','SUM']].to_csv('global_lineage_distribution_all.tsv',sep='\t')
	lineages = lineages[lineages['SUM'] >= 10].copy()
	lineages[['lineage1','lineage2','lineage3','lineage4','lineage5','lineage6','lineage7','SUM']].to_csv('global_lineage_distribution_filtred.tsv',sep='\t')



############
#DATA and Global Lineage Distribution Analysis
#
#############

def lineage_distribution_analysis():
	#What proportion of UK isolates (n=?) belong to lineage 1 and 3
	lineages = df.copy()
	lineages = lineages.replace({-1:0})
	UK_information = lineages[lineages['country'] == 'United Kingdom']
	number_UK = len(UK_information['strain'].values)
	number_UK_lineage13 = len(UK_information[(UK_information['lineage1'] == 1) | (UK_information['lineage3'] == 1)]['strain'].values)
	print("Proportion UK isolates belong lineage 1/3 is {}/{}: {} {}".format(number_UK_lineage13, number_UK, number_UK_lineage13/number_UK, cal(number_UK_lineage13, number_UK)))

	#Rest of analysis we are only looking at countries with 10 or more isolates 
	lineages['SUM'] = 1
	lineages = lineages.replace({-1:0}).groupby('country').sum().reset_index()
	lineages = lineages[lineages['SUM'] >= 10].copy()

	#Number of Lineage 1,2,3 in European isolates
	European = lineages[lineages['country'].isin(country_to_region['Europe'])]
	number_European = European['SUM'].sum()
	number_European_lineage123 = European['lineage1'].sum()+European['lineage2'].sum()+European['lineage3'].sum()

	print("Proportion European isolates lineage 123 is {}\t{}: {}".format(number_European_lineage123, number_European, 
		number_European_lineage123/number_European))


	#Number of Lineage 1,2,3 in North and South American isolates isolates
	NSAmerica = lineages[(lineages['country'].isin(country_to_region['North America'])) | (lineages['country'].isin(country_to_region['South America']))]
	number_NSAmerica = NSAmerica['SUM'].sum()
	number_NSAmerica_lineage123 = NSAmerica['lineage1'].sum()+NSAmerica['lineage2'].sum()+NSAmerica['lineage3'].sum()

	print("Proportion NSAmerica isolates lineage 123 is {}\t{}: {}".format(number_NSAmerica_lineage123, number_NSAmerica, 
		number_NSAmerica_lineage123/number_NSAmerica))

	#Number of Lineage 4 in Asia
	Asia = lineages[lineages['country'].isin(country_to_region['Asia'])]
	number_Asia = Asia['SUM'].sum()
	number_Asia_lineage4 = Asia['lineage4'].sum()

	print("Proportion Asian isolates lineage 4 is {}\t{}: {}".format(number_Asia_lineage4, number_Asia, 
		number_Asia_lineage4/number_Asia))

	#Number of Lineage 4 in Africa
	Africa = lineages[lineages['country'].isin(country_to_region['Africa'])]
	number_Africa = Africa['SUM'].sum()
	number_Africa_lineage4 = Africa['lineage4'].sum()

	print("Proportion Africa isolates lineage 4 is {}\t{}: {}".format(number_Africa_lineage4, number_Africa, 
		number_Africa_lineage4/number_Africa))



###########PERFORMING THE RESISTANCE DISTRIBUTION ANALYSIS ANALYSIS###############
############
#Phenotypic Resistant Distribution Analysis
#
#############
def generate_other_res(other_res):
	other_res = pd.DataFrame(list(other_res.items()), columns = ['Drug','Count'])
	other_res.to_csv('other_resistance.tsv',sep='\t')


def resistance_distribution_analysis():
	#How many  isolates and countires left?
	rd = df.copy()
	rd  = rd.replace({-1:0})

	#Only consider countries with more than 10 isolates WITH RESISTANT DATA and no "not provided" country
	rd = rd[(rd['NO_DATA'] != 1)&(rd['country'] != 'Not Provided')]
	rd =  rd.groupby('country').filter(lambda x: x['country'].count() >= 10)
	#Number of isolates and unique countries
	print("Number isolates with resistant data and from country with more than 10 isolates: {}".format(rd['strain'].count()))
	print("Number countries with more than 10 isolates with resistant data: {}".format(rd['country'].nunique()))

	#Number of isolates considered in resistance analysis
	print("Number of isolates considered {}".format(rd['strain'].nunique()))
	#Number of MDR isolates
	print("Number of MDR isolates: {}".format(rd[(rd['MDR'] == 1)]['strain'].nunique()))

	#Number of XDR isolates 
	print("Number of XDR isolates: {}".format(rd[(rd['XDR'] == 1)]['strain'].nunique()))

	#Number of isolates pan susceptible
	print("Number of S isolates: {}".format(rd[(rd['S'] == 1)]['strain'].nunique()))

	#Number of drug-resistant isolates 
	print("Number of DR isolates: {}".format(rd[(rd['S'] == 0)]['strain'].nunique()))

	#Number of countries with MDR?
	print("Number countries with MDR: {}".format(rd[rd['MDR'] == 1]['country'].nunique()))

	#Number of countries with MDR and XDR?
	print("Number of countreis with MDR and XDR: {}".format(rd[(rd['MDR'] == 1)&(rd['XDR'] == 1)]['country'].nunique()))

	#WHO_incidence_comparison()


	##Perform lineage MDR rate analysis
	print('LINEAGE1: {}\t{}'.format(rd[(rd['lineage1'] == 1)&(rd['MDR'] == 1)]['strain'].count(),rd[rd['lineage1'] == 1]['strain'].count()))
	print('LINEAGE2: {}\t{}'.format(rd[(rd['lineage2'] == 1)&(rd['MDR'] == 1)]['strain'].count(),rd[rd['lineage2'] == 1]['strain'].count()))
	print('LINEAGE3: {}\t{}'.format(rd[(rd['lineage3'] == 1)&(rd['MDR'] == 1)]['strain'].count(),rd[rd['lineage3'] == 1]['strain'].count()))
	print('LINEAGE4: {}\t{}'.format(rd[(rd['lineage4'] == 1)&(rd['MDR'] == 1)]['strain'].count(),rd[rd['lineage4'] == 1]['strain'].count()))
	##


	# test = rd.replace({-1:0}).groupby('country').sum()
	# test['SUM_W_DATA'] = test['SUM'] - test['NO_DATA']
	interest = ['S','MDR','XDR','INH_MONO', 'STR_MONO', 'OTHER_R','NO_DATA']
	# rd[interest].to_csv('global_resistance_distribution_all.tsv',sep='\t')
	# to_keep = test[interest].copy()
	# test = test[test['SUM_W_DATA'] >= 10].copy()
	rd[interest].to_csv('global_resistance_distribution_filtred.tsv',sep='\t')


	###########PERFORMING THE OHTER_R ANALYSIS###############
	other_res = {}
	def other_resistance(row):
		global count
		if(row['OTHER_R']):
			init = ''
			for drug in ['ISONIAZID','RIF', 'ETHAMBUTOL', 'PYRAZINAMIDE','STREPTOMYCIN','FQ','SLIS','ETH','CYS','PAS']:
				if(row[drug]):
					init = init + drug + '_'
			other_res[init] = other_res.get(init,0) + 1


	#Of the total number of DR isolates how many had other_r
	print("Number DR isolates {}".format(rd[rd['S'] != 1]['strain'].count()))
	print("Number Susceptible isolates {}".format(rd[rd['S'] == 1]['strain'].count()))
	print("Number of Other_R isolates {}".format(rd[rd['OTHER_R'] == 1]['strain'].count()))


	#Find proportion of lineage with this other resistance
	secret = rd.copy()
	sub = secret[secret['OTHER_R'] == 1].sum()[['lineage1','lineage2','lineage3','lineage4']].copy()
	tot = secret.sum()[['lineage1','lineage2','lineage3','lineage4']].copy()
	print(sub)
	print(tot)
	print(sub/tot)

	#Among top five countries with data, country with max number of other resistance
	secret = rd.copy()
	top_five_countries = [i[0] for i in secret.groupby('country').country.value_counts().nlargest(5).index.tolist()]
	print(top_five_countries)
	secret = secret[secret['country'].isin(top_five_countries)]
	# filter(lambda x: x['country'].count() >= 10)
	# secret = secret 
	sub = secret[secret['OTHER_R'] == 1].groupby('country').count()['strain'].copy()
	tot = secret.groupby('country').count()['strain'].copy()
	print(sub)
	print(tot)
	print(((sub/tot).nlargest(4)))


	rd.apply(lambda row: other_resistance(row),axis=1)
	max_category = max(other_res.items(), key=operator.itemgetter(1))[0]
	second_max_category = max([i for i in other_res.items() if max_category != i[0]], key=operator.itemgetter(1))[0]
	print("Most other {} at {}/{} {}".format(max_category, other_res[max_category],rd[rd['OTHER_R'] == 1]['strain'].count(),other_res[max_category]/rd[rd['OTHER_R'] == 1]['strain'].count()))
	print("Second most  other {} at {}/{} {}".format(second_max_category, other_res[second_max_category],rd[rd['OTHER_R'] == 1]['strain'].count(),other_res[second_max_category]/rd[rd['OTHER_R'] == 1]['strain'].count()))

	number_flq = 0
	for i in other_res:
		if 'FQ' in i:
			number_flq += other_res[i]

	print("FLQ NUMBER {}/{} {}".format(number_flq,rd[rd['OTHER_R'] == 1]['strain'].count(),number_flq/rd[rd['OTHER_R'] == 1]['strain'].count()))
	#generate_other_res(other_res)

# resistance_distribution_analysis()

def WHO_incidence_comparison():
	print("WHO INCIDENCE COMPARISON ANALYSIS")
	#How many  isolates and countires left?
	rd = df.copy()
	rd  = rd.replace({-1:0})

	#Only consider countries with more than 10 isolates WITH RESISTANT DATA and no "not provided" country
	rd = rd[(rd['NO_DATA'] != 1)&(rd['country'] != 'Not Provided')]
	rd =  rd.groupby('country').filter(lambda x: x['country'].count() >= 10)
	greater = 0
	equal = 0
	less = 0


	for country in WHO_incidence:
		start_WHO,end_WHO = WHO_incidence[country][1], WHO_incidence[country][2]
		country_df = rd[rd['country'] == country].copy()
		number = country_df['strain'].count()
		number_MDR = country_df[country_df['MDR'] == 1]['strain'].count()
		prop_MDR = number_MDR/number 
		start_exp, end_exp = sp.proportion_confint(number_MDR, number)
		if(prop_MDR != 0.0):	
			if(end_exp < start_WHO):
				less += 1
				print("{}\t{}".format(country, 'LESS'))
			elif(start_exp > end_WHO):
				greater += 1
				print("{}\t{}".format(country, 'GREATER'))
			else:
				equal += 1
				print("{}\t{}".format(country, 'EQUAL'))

	print("GREATER: {}".format(greater))
	print("EQUAL: {}".format(equal))
	print("LESS: {}".format(less))
	print("TOTAL: {}".format(greater+equal+less))

# resistance_distribution_analysis()

# ###########PERFORMING THE MDR RECENT AQUISITION###############
def relative_order_of_phenotypic_resistance():
	drugs_present = [i for i in list(set(name_with_date['drug'].values)) if i not in ['PAS','NICOTINAMIDE']]
	for drug in drugs_present:
		IQR = name_with_date[name_with_date['drug'] == drug]['date'].quantile([0.25,0.5,0.75])
		print("DRUG: {} STATS: {}".format(drug, IQR))

	groups = []
	for drug_one in drugs_present:
		drug_one_added = False
		for drug_two in drugs_present:
			if(drug_one != drug_two):
				one = name_with_date[(name_with_date['drug'] == drug_one)]['date'].copy()
				two = name_with_date[(name_with_date['drug'] == drug_two)]['date'].copy()
				print("{}\t{}\t{}".format(drug_one, drug_two,sc.ranksums(one,two).pvalue))
				if(sc.ranksums(one,two).pvalue >= 0.01):
					drug_one_added = True
					add_group = False
					for group in groups:
						if(drug_one in group and drug_two not in group):
							group.append(drug_two)
							add_group = False
						elif(drug_two in group and drug_one not in group):
							group.append(drug_one)
							add_group = False
						elif(drug_one in group and drug_two in group):
							add_group = False
							break
						else:
							add_group = True
					if(add_group or not len(groups)):
						groups.append([drug_one, drug_two])
		if not drug_one_added:
			groups.append([drug_one])

	print(groups)


# relative_order_of_phenotypic_resistance()

def median_MRSCA_correlate_dateintro():
	x_value = []
	y_value = []

	for drug in date_drug_introduced:
		x_value.append(date_drug_introduced[drug])
		y_value.append(name_with_date[name_with_date['drug'] == drug]['date'].median())
	print(x_value)
	print(y_value)
	result = sm.OLS(y_value, sm.add_constant(x_value)).fit()
	print(result.summary())
	p = result.params 
	x = np.array(x_value)
	plt.scatter(x_value, y_value)
	print(result.rsquared)
	plt.plot(x_value, p[0] + p[1]*x)
	plt.show()

#median_MRSCA_correlate_dateintro()
	


def recent_aquisition():
	dates = df.copy()
	dates = dates.replace({-1:0}).copy()
	to_keep = dates.groupby('country').sum().copy()
	del dates['date']
	dates= pd.merge(dates, name_with_date, on='strain',how='outer')
	test = dates[(dates['MDR'] == 1) & (dates['date'] <= 5)].copy()
	test['Count'] = 1
	sub_test = test.groupby(['country','drug']).sum().reset_index()[['country','drug','Count']].copy()
	sub_result = sub_test.merge(to_keep,on='country')[['country','drug','Count','MDR']].copy()
	sub_result['percent of MDR isolates'] = sub_result['Count']/sub_result['MDR']
	sub_result['lower_bound'] = sub_result.apply(lambda row: confidence_interval_lower(row.Count,row.MDR),axis=1)
	sub_result['upper_bound'] = sub_result.apply(lambda row: confidence_interval_upper(row.Count, row.MDR),axis=1)
	sub_result = sub_result.rename(columns={'count': 'number of MDR isolates', 'MDR': 'number of MDR isolates in country'})
	#sub_result.to_csv('recent_aquistion.tsv',sep='\t')
	#print(sub_result)

	#How many countries are here?
	print("Number of countries {}".format(sub_result['country'].nunique()))

	#How many countries have measurable for PZA and EMB
	PEMB = sub_result[(sub_result['drug'] == 'PZA') | (sub_result['drug'] == 'EMB')].copy()
	PEMB_confidence_measurable = PEMB[PEMB['lower_bound'].astype(float) >= 0.01]
	country_PEMB_count = PEMB_confidence_measurable.replace({'PZA': 1, 'EMB':1}).groupby('country').sum().reset_index()
	countries_with_both = country_PEMB_count[country_PEMB_count['drug'] == 2]['country'].values
	PEMB = PEMB_confidence_measurable[PEMB_confidence_measurable['country'].isin(countries_with_both)]
	print("Countries with measurable for PZA and EMB {}".format(countries_with_both))
	max_PEMB = PEMB.loc[PEMB['percent of MDR isolates'].idxmax()]
	print("MAX measurable PZA and EMB is {}\t{}({}-{})".format(max_PEMB['country'], max_PEMB['percent of MDR isolates'], max_PEMB['lower_bound'], max_PEMB['upper_bound']))
	min_PEMB = PEMB.loc[PEMB['percent of MDR isolates'].idxmin()]
	print("MIN measurable PZA and EMB is {}\t{}({}-{})".format(min_PEMB['country'], min_PEMB['percent of MDR isolates'], min_PEMB['lower_bound'], min_PEMB['upper_bound']))


	#How many countries have measurable for FLQ and SLI
	FQSL = sub_result[(sub_result['drug'] == 'FLQ') | (sub_result['drug'] == 'SLIs')].copy()
	FQSL_confidence_measurable = FQSL[FQSL['lower_bound'].astype(float) >= 0.01]
	country_FQSL_count = FQSL_confidence_measurable.replace({'FLQ': 1, 'SLIs':1}).groupby('country').sum().reset_index()
	countries_with_both = country_FQSL_count[country_FQSL_count['drug'] == 2]['country'].values
	FQSL = FQSL_confidence_measurable[FQSL_confidence_measurable['country'].isin(countries_with_both)]
	print("Countries with measurable for FQL and SLI {}".format(countries_with_both))
	max_FQSL_FQ = FQSL[FQSL['drug'] == 'FLQ'].loc[FQSL[FQSL['drug'] == 'FLQ']['percent of MDR isolates'].idxmax()]
	print("MAX measurable FQL is {}\t{}({}-{})".format(max_FQSL_FQ['country'], max_FQSL_FQ['percent of MDR isolates'], max_FQSL_FQ['lower_bound'], max_FQSL_FQ['upper_bound']))
	min_FQSL_FQ = FQSL[FQSL['drug'] == 'FLQ'].loc[FQSL[FQSL['drug'] == 'FLQ']['percent of MDR isolates'].idxmin()]
	print("MIN measurable FQL is {}\t{}({}-{})".format(min_FQSL_FQ['country'], min_FQSL_FQ['percent of MDR isolates'], min_FQSL_FQ['lower_bound'], min_FQSL_FQ['upper_bound']))
	max_FQSL_SL = FQSL[FQSL['drug'] == 'SLIs'].loc[FQSL[FQSL['drug'] == 'SLIs']['percent of MDR isolates'].idxmax()]
	print("MAX measurable SLI is {}\t{}({}-{})".format(max_FQSL_SL['country'], max_FQSL_SL['percent of MDR isolates'], max_FQSL_SL['lower_bound'], max_FQSL_SL['upper_bound']))
	min_FQSL_SL = FQSL[FQSL['drug'] == 'SLIs'].loc[FQSL[FQSL['drug'] == 'SLIs']['percent of MDR isolates'].idxmin()]
	print("MIN measurable SLI is {}\t{}({}-{})".format(min_FQSL_SL['country'], min_FQSL_SL['percent of MDR isolates'], min_FQSL_SL['lower_bound'], min_FQSL_SL['upper_bound']))


	#median MRSCA age for FLQ or SLI resistance acquisition among MDR isolates
	dates = dates[(dates['MDR'] == 1)].copy()
	# MFQSL = dates[(dates['FQ'] == 1) | (dates['SLIS'] == 1)]
	#No! the above is wrong we have to use the dates['drug'] since that is teh drug that was dated! other way could give us wrong dates
	MFQSL = dates[(dates['drug'] == 'FLQ') | (dates['drug'] == 'SLIs')]
	print("MEDIAN MRSCA age for FLQ or SLI res aquisition among MDR isolates {}".format(MFQSL['date'].quantile([0.25,0.5,0.75])))

# recent_aquisition()

#####OOLD RESISTANCE STUFFFFF#######

def old_resistance():
	dates = pd.read_csv('strain_name_with_date_unprocessed',sep='\t')[['country','strain','drug','date']].copy()
	# dates['drug'] = dates['drug'].replace({'PYRAZINAMIDE':'PZA','ISONIAZID':'INH','FLUOROQUINOLONES':'FLQ', 'STREPTOMYCIN':'STR',
	# 	'RIFAMPICIN/RIFABUTIN':'RIF','AMINOGLYCOSIDE':'SLIs','ETHAMBUTOL':'EMB'})
	dates['drug']=dates['drug'].replace({'MOXIFLOXACIN':'FLQ','CIPROFLOXACIN':'FLQ','OFLOXACIN':'FLQ'})
	dates['drug']=dates['drug'].replace({'KANAMYCIN':'SLIs','AMIKACIN':'SLIs','CAPREOMYCIN':'SLIs'})
	dates['drug']=dates['drug'].replace({'RIFAMPICIN':'RIF','RIFABUTIN':'RIF','RIFAMPICIN/RIFABUTIN':'RIF'})
	dates['drug']=dates['drug'].replace({'ETHIONAMIDE':'ETH','PROTHIONAMIDE':'ETH'})
	dates['drug']=dates['drug'].replace({'CYCLOSERINE':'CYS'})
	dates['drug']=dates['drug'].replace({'PARA-AMINOSALICYLIC_ACID':'PAS'})
	dates['drug']=dates['drug'].replace({'PYRAZINAMIDE':'PZA'})
	dates['drug']=dates['drug'].replace({'STREPTOMYCIN':'STR'})
	dates['drug']=dates['drug'].replace({'ISONIAZID':'INH'})
	oldRifdenom = dates[(dates['drug'] == 'RIF')]['date'].count()
	oldRifnum = dates[(dates['drug'] == 'RIF')&(dates['date'] >= 20)]['date'].count()
	oldINHdenom = dates[(dates['drug'] == 'INH')]['date'].count()
	oldINHnum = dates[(dates['drug'] == 'INH')&(dates['date'] >= 20)]['date'].count()
	oldFLQdenom = dates[(dates['drug'] == 'FLQ')]['date'].count()
	oldFLQnum = dates[(dates['drug'] == 'FLQ')&(dates['date'] >= 20)]['date'].count()
	oldRIF = oldRifnum/oldRifdenom
	oldINH = oldINHnum/oldINHdenom
	oldFLQ = oldFLQnum/oldFLQdenom
	print("Percentage of old RIF {} (95% CI {}-{})".format(oldRIF*100, confidence_interval_lower(oldRifnum, oldRifdenom)*100, confidence_interval_upper(oldRifnum, oldRifdenom)*100))
	print("Fraction for old RIF {}/{}".format(oldRifnum, oldRifdenom))
	print("Percentage of old INH {} (95% CI {}-{})".format(oldINH*100, confidence_interval_lower(oldINHnum, oldINHdenom)*100, confidence_interval_upper(oldINHnum, oldINHdenom)*100))
	print("Fraction for old INH {}/{}".format(oldINHnum, oldINHdenom))
	print("Percentage of old FLQ {} (95% CI {}-{})".format(oldFLQ*100, confidence_interval_lower(oldFLQnum, oldFLQdenom)*100, confidence_interval_upper(oldFLQnum, oldFLQdenom)*100))
	print("Fraction for old FLQ {}/{}".format(oldFLQnum, oldFLQdenom))

	###how many countries for each had old resistance
	def countries(drug):
		country_counts = dates[dates['drug'] == drug].groupby('country').count().reset_index()
		country_counts = dict(zip(country_counts.country, country_counts.strain))
		#oldrugCountrydenom = len(set([i for i in dates[(dates['drug'] == drug)]['country'].values if country_counts[i] > 10]))
		#oldrugnum = len(set([i for i in dates[(dates['drug'] == drug)&(dates['date'] >= 20)]['country'].values if country_counts[i] > 10]))
		oldrugCountrydenom = len(set([i for i in dates[(dates['drug'] == drug)]['country'].values]))
		oldrugnum = len(set([i for i in dates[(dates['drug'] == drug)&(dates['date'] >= 20)]['country'].values]))
		return oldrugCountrydenom, oldrugnum

	oldRifCountrydenom, oldRifnum = countries('RIF')
	oldINHCountrydenom, oldINHnum = countries('INH')
	oldFLQCountrydenom, oldFLQnum = countries('FLQ')
	print("Percentage of old RIF {}, {}/{})".format(oldRifnum, oldRifnum, oldRifCountrydenom))
	print("Percentage of old INH {}, {}/{})".format(oldINHnum, oldINHnum, oldINHCountrydenom))
	print("Percentage of old FLQ {}, {}/{})".format(oldFLQnum, oldFLQnum, oldFLQCountrydenom))

def geographic_distribution_MRSCA_comparison():
	dates = df.copy()
	dates = dates.replace({-1:0}).copy()
	to_keep = dates.groupby('country').sum().copy()
	del dates['date']
	#name_with_date = pd.read_csv('strain_name_with_date.tsv',sep='\t')[['strain','drug','date']].copy()
	dates= pd.merge(dates, name_with_date, on='strain',how='inner')
	four_drug_interest = ['INH','RIF','SLIs','FLQ']
	five_countries_most_resistance = dates[(dates['S'] != 1)&(dates['country'] != 'Not Provided')&(dates['NO_DATA'] != 1)].groupby('country').count().nlargest(5, 'strain').reset_index()['country'].values
	print(five_countries_most_resistance)

	for drug in ['INH', 'RIF','SLIs','PZA','ETH','FLQ']:
		if(drug == 'FLQ'):
			five_countries_most_resistance = [i for i in five_countries_most_resistance if i != 'Russia']
		groups = []
		for country_one in five_countries_most_resistance:
			country_one_added = False
			for country_two in five_countries_most_resistance:
				if(country_one != country_two):
					one = dates[(dates['drug'] == drug)&(dates['country'] == country_one)]['date']
					two = dates[(dates['drug'] == drug)&(dates['country'] == country_two)]['date']
					print("{}\t{}\t{}".format(country_one, country_two,sc.ranksums(one,two).pvalue))
					if(sc.ranksums(one,two).pvalue >= 0.01):
						country_one_added = True
						add_group = False
						for group in groups:
							if(country_one in group and country_two not in group):
								group.append(country_two)
								add_group = False
							elif(country_two in group and country_one not in group):
								group.append(country_one)
								add_group = False
							elif(country_one in group and country_two in group):
								add_group = False
								break
							else:
								add_group = True
						if(add_group or not len(groups)):
							groups.append([country_one, country_two])
			if not country_one_added:
				groups.append([country_one])
		print("GROUP FOR {} \t {}".format(drug, groups))
		for group in groups:
			print("IQR for {} is {}".format(group, dates[(dates['drug'] == drug)&(dates['country'].isin(group))]['date'].quantile([0.25,0.5,0.75])))

	for drug in ['RIF','SLIs','FLQ']:
		print("IQR for United Kingdom {} is {}".format(drug, dates[(dates['drug'] == drug)&(dates['country'].isin(['United Kingdom']))]['date'].quantile([0.25,0.5,0.75])))

# geographic_distribution_MRSCA_comparison()

def MRSCA_ratio():
	test = df.copy()
	test = test.replace({-1:0}).copy()
	test = pd.merge(test, name_with_date, on='strain',how='inner')
	test['SUM'] = 1
	final = {}

	####FOR POOLED ACROSS COUNTRY########

	#How many dated strains in each country?
	country_total = test.groupby('country').sum()
	country_total = country_total[country_total['SUM'] > 10]
	country_total = country_total['SUM']

	#divide for answer
	result = '('+test.round(3).groupby('country')['date_y'].nunique().astype('str')+'/'+country_total.astype('str')+')'+ (1-(test.round(3).groupby('country')['date_y'].nunique()/country_total).round(2)).round(2).astype('str')
	final['pooled'] = result.to_dict()


	####FOR SPECIFIC DRUGS#####

	for drug in ['INH','RIF','PZA','FLQ','STR','SLIs','EMB']:
		#How many dated strains in each country?
		country_total = test[test['drug'] == drug].groupby('country').sum()
		country_total = country_total[country_total['SUM'] > 10]
		country_total = country_total['SUM']
		#divide for answer
		result = '('+test[test['drug'] == drug].round(3).groupby('country')['date_y'].nunique().astype('str')+'/'+country_total.astype('str')+')'+ (1-(test[test['drug'] == drug].round(3).groupby('country')['date_y'].nunique()/country_total).round(2)).round(2).astype('str')
		final[drug] = result

	resultat = pd.DataFrame.from_dict(final)
	#resultat.to_csv('MRSCA_ratios.tsv',sep='\t')
	print(resultat)

#####DISTRIBUTION OF RESISTANCE MUTATIONS SECTION#######


def calculate_resistant_mutations_stats():
	strain_info = df.copy()

	resistance_mutation = pd.read_csv('results_modified_unknown',sep='\t')
	combined = pd.merge(strain_info, resistance_mutation, on='strain', how='inner')

	##Since we are considering variation over geography we need to filter out isolates without a geography
	combined = combined[combined['country'] != 'Not Provided']

	del combined['Unnamed: 0']

	#Create a lineage column to make lineage pooling easier
	def lineage(row):
		if(row['lineage1']):
			return 'lineage1'
		elif(row['lineage2']):
			return 'lineage2'
		elif(row['lineage3']):
			return 'lineage3'
		elif(row['lineage4']):
			return 'lineage4'

	combined['lineage'] = combined.apply(lambda row: lineage(row), axis=1)

	#Step two now group based on country
	combined = combined.replace({-1:0})
	combined['SUM'] = 1
	sub = combined.groupby('country').sum().reset_index()


	#Step 3 now divide each mutation column by the total number of resistant isolates
	mutations = [i for i in sub.columns if '_' in i and i not in ['AMOXICILLIN_CLAVULANATE', 'PARA-AMINOSALICYLIC_ACID', 'NO_DATA'] and 'unknown' not in i and 'MONO' not in i and 'OTHER' not in i]
	for mutation in mutations:
		drug = drug_mapping[mutation.split('_')[0]]
		filtered = combined[combined[drug] == 1].groupby('country').sum().reset_index()
		filtered = filtered[filtered['SUM'] > 10]
		sub[mutation] = filtered[mutation]/filtered['SUM']


	#Step 4 now create a variance per mutation row for each mutation
	mutation_to_variance = {}
	for mutation in mutations:
		ugh = [i for i in list(sub[mutation].values) if str(i) != 'nan' and str(i) != 'inf']
		mutation_to_variance[mutation] = np.std(ugh)
		if(('rpoB' in mutation and 'I491F' in mutation)):
			ugh = [i for i in list(sub[mutation].values) if str(i) != 'nan' and str(i) != 'inf']
			print("stats for rpoB I491F: {}\t{}\t{}-{}".format(mutation, np.std(ugh), min(ugh), max(ugh)))


	#Step 5 output 10 variants with largest variance
	result = dict(sorted(mutation_to_variance.items(), key=operator.itemgetter(1), reverse=True)[:10])
	for mutation in result.keys():
		ugh = [i for i in list(sub[mutation].values) if str(i) != 'nan' and str(i) != 'inf']
		print('FOR COUNTRY: {}\t{}\t{}-{}'.format(mutation, np.std(ugh), min(ugh), max(ugh)))


	#Step 6 find which of original 267 resistant mutations has higher variance than INH_SNP_P_1673425_CT.15_fabG1.inhA 
	original_snps = [i.rstrip().replace('\t','_') for i in open('snps','r').readlines()]
	CT_15_std = mutation_to_variance['INH_SNP_P_1673425_CT.15_fabG1.inhA']
	ugh = [i for i in list(sub['INH_SNP_P_1673425_CT.15_fabG1.inhA'].values) if str(i) != 'nan' and str(i) != 'inf']
	print("fabG1 info: {}\t{}-{}".format(np.std(ugh), np.min(ugh), np.max(ugh)))
	number_greater = 0
	for i in original_snps:
		if(mutation_to_variance[i] > CT_15_std):
			number_greater += 1
			ugh = [i for i in list(sub[i].values) if str(i) != 'nan' and str(i) != 'inf']
			print('GREATER: {}\t{}%\t{}%-{}%'.format(i, round(np.std(ugh)*100,2), round(min(ugh)*100,2), round(max(ugh)*100,2)))
	print("NUMBER GREATER THAN CT 15 {}/{}".format(number_greater, len(original_snps)))


	#LINEAGE Step two now group based on lineage
	combined = combined.replace({-1:0})
	combined['SUM'] = 1
	sub = combined.groupby('lineage').sum().reset_index()


	#Step 3 now divide each mutation column by the total number of resistant isolates
	mutations = [i for i in sub.columns if '_' in i and i not in ['AMOXICILLIN_CLAVULANATE', 'PARA-AMINOSALICYLIC_ACID', 'NO_DATA'] and 'unknown' not in i and 'MONO' not in i and 'OTHER' not in i]
	for mutation in mutations:
		drug = drug_mapping[mutation.split('_')[0]]
		filtered = combined[combined[drug] == 1].groupby('lineage').sum().reset_index()
		filtered = filtered[filtered['SUM'] > 10]
		sub[mutation] = filtered[mutation]/filtered['SUM']

	#Step 4 now create a variance per mutation row for each mutation
	mutation_to_variance = {}
	for mutation in mutations:
		ugh = [i for i in list(sub[mutation].values) if str(i) != 'nan' and str(i) != 'inf']
		mutation_to_variance[mutation] = np.std(ugh)


	#Step 5 output 10 variants with largest variance
	result = dict(sorted(mutation_to_variance.items(), key=operator.itemgetter(1), reverse=True)[:10])
	for mutation in result.keys():
		ugh = [i for i in list(sub[mutation].values) if str(i) != 'nan' and str(i) != 'inf']
		print('FOR LINEAGE: {}\t{}\t{}-{}'.format(mutation, np.std(ugh), min(ugh), max(ugh)))


	#Step 6 find which of original 267 resistant mutations has higher variance than INH_SNP_P_1673425_CT.15_fabG1.inhA 
	original_snps = [i.rstrip().replace('\t','_') for i in open('snps','r').readlines()]
	CT_15_std = mutation_to_variance['INH_SNP_P_1673425_CT.15_fabG1.inhA']
	ugh = [i for i in list(sub['INH_SNP_P_1673425_CT.15_fabG1.inhA'].values) if str(i) != 'nan' and str(i) != 'inf']
	print("fabG1 info: {}\t{}-{}".format(np.std(ugh), np.min(ugh), np.max(ugh)))
	number_greater = 0
	for i in original_snps:
		if(mutation_to_variance[i] > CT_15_std):
			number_greater += 1
			ugh = [i for i in list(sub[i].values) if str(i) != 'nan' and str(i) != 'inf']
			print('GREATER: {}\t{}%\t{}%-{}%'.format(i, round(np.std(ugh)*100,2), round(min(ugh)*100,2), round(max(ugh)*100,2)))
	print("NUMBER GREATER THAN CT 15 {}/{}".format(number_greater, len(original_snps)))

calculate_resistant_mutations_stats()


# ###########DOING CALCULATIONS FOR PAPER ABOUT MRSCA#######################

# """
# median MRSCA for RIF resistance acquisition among INH resistant isolates is X (HPD Y-Z), 
# and median MRSCA for FLQ or SLI resistance acquisition among MDR isolates was A (HPD B-C), 
# with a difference of E years
# """
# import pymc3

# test = df.copy()
# test = test.replace({-1:0}).copy()
# # test.apply(lambda row: other_resistance(row),axis=1)
# del test['date']
# test = pd.merge(test, name_with_date, on='strain',how='outer')

# median_MRSCA_RIF_among_INH = test[(test['ISONIAZID'] == 1)&(test['RIF'] == 1)]['date'].dropna().median()
# print(median_MRSCA_RIF_among_INH)
# print(pymc3.stats.hpd(test[(test['ISONIAZID'] == 1)&(test['RIF'] == 1)]['date'].dropna()))

# median_MRSCA_FLQorSLI_among_MDR = test[(test['MDR'] == 1)&((test['FQ'] == 1) | (test['SLIS'] == 1))]['date'].dropna().median()
# print(median_MRSCA_FLQorSLI_among_MDR)
# print(pymc3.stats.hpd(test[(test['MDR'] == 1)&((test['FQ'] == 1) | (test['SLIS'] == 1))]['date'].dropna()))

# print(median_MRSCA_RIF_among_INH-median_MRSCA_FLQorSLI_among_MDR)






