"""
Written by Yasha Ektefaie
August 2018

python3 results_pre_rm.py 

Run in a directory with all the data and depending on which method you call outputs the statistics for the results
section of the manuscript all before the resistance mutation part

"""

import pandas as pd 
import numpy as np  
import statsmodels.stats.proportion as sp
import operator
import scipy.stats as sc
import itertools

#Begin by taking the strain_info file and annotating it in the support function/def section
df = pd.read_csv('strain_info.tsv',sep='\t')
strains_we_have = [i.rstrip() for i in open('final_all_strains','r').readlines()]
df = df[df['strain'].isin(strains_we_have)]

##########Support Functions and Definitions#########


#Calculating confidence interval functions
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



######################################################

"""
LINEAGE ANALYSIS

This next section of functions have the code used to produce the results
for the "Data and global lineage distribution" section of the manuscript

"""

#Function to generate how many of each lineage were there for each country filtered and not filtered for countries with less than 10 isoaltes
def generate_lineage_data():
	lineages = df.copy()
	lineages['SUM'] = 1
	lineages = lineages.replace({-1:0}).groupby('country').sum()
	lineages[['lineage1','lineage2','lineage3','lineage4','lineage5','lineage6','lineage7','SUM']].to_csv('global_lineage_distribution_all.tsv',sep='\t')
	lineages = lineages[lineages['SUM'] >= 10].copy()
	lineages[['lineage1','lineage2','lineage3','lineage4','lineage5','lineage6','lineage7','SUM']].to_csv('global_lineage_distribution_filtred.tsv',sep='\t')


#Function that answers specific questions about the lineage distribution asked in manuscript
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

"""
RESISTANCE DISTRIBUTION ANALYSIS

This next section of functions have the code used to produce the results
for the "Phenotypic resistance distribution" section of the manuscript

"""

#Function that produces the category of "other resistant" strains and the number of isolates that fall into each cateogry
def generate_other_res(other_res):
	other_res = pd.DataFrame(list(other_res.items()), columns = ['Drug','Count'])
	other_res.to_csv('other_resistance.tsv',sep='\t')


#Function that answers specific questions about the resistant distribution asked in manuscript
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

	#Number of countries with MDR?
	print("Number countries with MDR: {}".format(rd[rd['MDR'] == 1]['country'].nunique()))

	#Number of countries with MDR and XDR?
	print("Number of countreis with MDR and XDR: {}".format(rd[(rd['MDR'] == 1)&(rd['XDR'] == 1)]['country'].nunique()))


	#######For this next part of the analysis we don't care about excluding countries with less than 10 isolates######

	test = df.copy()
	test['SUM'] = 1
	test = test.replace({-1:0})
	##Perform lineage MDR rate analysis
	print('LINEAGE1: {}\t{}'.format(test[test['lineage1'] == 1].sum()['MDR'],test[test['lineage1'] == 1].sum()['SUM']))
	print('LINEAGE2: {}\t{}'.format(test[test['lineage2'] == 1].sum()['MDR'],test[test['lineage2'] == 1].sum()['SUM']))
	print('LINEAGE3: {}\t{}'.format(test[test['lineage3'] == 1].sum()['MDR'],test[test['lineage3'] == 1].sum()['SUM']))
	print('LINEAGE4: {}\t{}'.format(test[test['lineage4'] == 1].sum()['MDR'],test[test['lineage4'] == 1].sum()['SUM']))
	##


	test = test.replace({-1:0}).groupby('country').sum()
	test['SUM_W_DATA'] = test['SUM'] - test['NO_DATA']
	interest = ['S','MDR','XDR','INH_MONO', 'STR_MONO', 'OTHER_R','NO_DATA', 'SUM', 'SUM_W_DATA']
	test[interest].to_csv('global_resistance_distribution_all.tsv',sep='\t')
	to_keep = test[interest].copy()
	test = test[test['SUM_W_DATA'] >= 10].copy()
	test[interest].to_csv('global_resistance_distribution_filtred.tsv',sep='\t')


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

	test = df.copy()
	test = test.replace({-1:0}).copy()

	#Of the total number of DR isolates how many had other_r
	print("Number DR isolates {}".format(rd[rd['S'] != 1]['strain'].count()))
	print("Number Susceptible isolates {}".format(rd[rd['S'] == 1]['strain'].count()))
	print("Number of Other_R isolates {}".format(test[test['OTHER_R'] == 1]['strain'].count()))


	#Find proportion of lineage with this other resistance
	secret = test.copy()
	sub = secret[secret['OTHER_R'] == 1].sum()[['lineage1','lineage2','lineage3','lineage4']].copy()
	tot = secret.sum()[['lineage1','lineage2','lineage3','lineage4']].copy()
	print(sub/tot)

	#Find country with max number of other resistance
	secret = test.copy()
	sub = secret[secret['OTHER_R'] == 1].groupby('country').count()['strain'].copy()
	tot = secret.groupby('country').count()['strain'].copy()
	print((sub.nlargest(4)/tot[sub.nlargest(4).index]))

	#What country had the most and second most other_r resistant isolates
	test.apply(lambda row: other_resistance(row),axis=1)
	max_category = max(other_res.items(), key=operator.itemgetter(1))[0]
	second_max_category = max([i for i in other_res.items() if max_category != i[0]], key=operator.itemgetter(1))[0]
	print("Most other {} at {}/{} {}".format(max_category, other_res[max_category],test[test['OTHER_R'] == 1]['strain'].count(),other_res[max_category]/test[test['OTHER_R'] == 1]['strain'].count()))
	print("Second most  other {} at {}/{} {}".format(second_max_category, other_res[second_max_category],test[test['OTHER_R'] == 1]['strain'].count(),other_res[second_max_category]/test[test['OTHER_R'] == 1]['strain'].count()))

	number_flq = 0
	for i in other_res:
		if 'FQ' in i:
			number_flq += other_res[i]

	print("FLQ NUMBER {}/{} {}".format(number_flq,test[test['OTHER_R'] == 1]['strain'].count(),number_flq/test[test['OTHER_R'] == 1]['strain'].count()))
	generate_other_res(other_res)

"""
DATING ANALYSIS

This next section of functions have the code used to produce the results
for the "Molecular dating of Resistance Acquisition" section of the manuscript

"""

#Calculations for recent aqusition stats in paragraph that starts with "We assessed the frequency of recent resistance amplification" 
def recent_aquisition():
	dates = df.copy()
	dates = dates.replace({-1:0}).copy()
	to_keep = dates.groupby('country').sum().copy()
	del dates['date']
	name_with_date = pd.read_csv('strain_name_with_date.tsv',sep='\t')[['strain','drug','date']].copy()
	dates= pd.merge(dates, name_with_date, on='strain',how='outer')
	test = dates[(dates['MDR'] == 1) & (dates['date'] <= 5)].copy()
	test['Count'] = 1
	sub_test = test.groupby(['country','drug']).sum().reset_index()[['country','drug','Count']].copy()
	sub_result = sub_test.merge(to_keep,on='country')[['country','drug','Count','MDR']].copy()
	sub_result['percent of MDR isolates'] = sub_result['Count']/sub_result['MDR']
	sub_result['lower_bound'] = sub_result.apply(lambda row: confidence_interval_lower(row.Count,row.MDR),axis=1)
	sub_result['upper_bound'] = sub_result.apply(lambda row: confidence_interval_upper(row.Count, row.MDR),axis=1)
	sub_result = sub_result.rename(columns={'count': 'number of MDR isolates', 'MDR': 'number of MDR isolates in country'})

	####Uncomment this if you want to produce recent_aquisition distribution tsv to use in plotting
	#sub_result.to_csv('recent_aquistion.tsv',sep='\t')

	#How many countries are here?
	print("Number of countries {}".format(sub_result['country'].nunique()))

	#How many countries have measurable for PZA and EMB
	PEMB = sub_result[(sub_result['drug'] == 'PYRAZINAMIDE') | (sub_result['drug'] == 'ETHAMBUTOL')].copy()
	PEMB_confidence_measurable = PEMB[PEMB['lower_bound'].astype(float) >= 0.01]
	country_PEMB_count = PEMB_confidence_measurable.replace({'PYRAZINAMIDE': 1, 'ETHAMBUTOL':1}).groupby('country').sum().reset_index()
	countries_with_both = country_PEMB_count[country_PEMB_count['drug'] == 2]['country'].values
	PEMB = PEMB_confidence_measurable[PEMB_confidence_measurable['country'].isin(countries_with_both)]
	print("Countries with measurable for PZA and EMB {}".format(countries_with_both))
	max_PEMB = PEMB.loc[PEMB['percent of MDR isolates'].idxmax()]
	print("MAX measurable PZA and EMB is {}\t{}({}-{})".format(max_PEMB['country'], max_PEMB['percent of MDR isolates'], max_PEMB['lower_bound'], max_PEMB['upper_bound']))
	min_PEMB = PEMB.loc[PEMB['percent of MDR isolates'].idxmin()]
	print("MIN measurable PZA and EMB is {}\t{}({}-{})".format(min_PEMB['country'], min_PEMB['percent of MDR isolates'], min_PEMB['lower_bound'], min_PEMB['upper_bound']))


	#How many countries have measurable for FLQ and SLI
	FQSL = sub_result[(sub_result['drug'] == 'FLUOROQUINOLONES') | (sub_result['drug'] == 'AMINOGLYCOSIDE')].copy()
	FQSL_confidence_measurable = FQSL[FQSL['lower_bound'].astype(float) >= 0.01]
	country_FQSL_count = FQSL_confidence_measurable.replace({'FLUOROQUINOLONES': 1, 'AMINOGLYCOSIDE':1}).groupby('country').sum().reset_index()
	countries_with_both = country_FQSL_count[country_FQSL_count['drug'] == 2]['country'].values
	FQSL = FQSL_confidence_measurable[FQSL_confidence_measurable['country'].isin(countries_with_both)]
	print("Countries with measurable for FQL and SLI {}".format(countries_with_both))
	max_FQSL_FQ = FQSL[FQSL['drug'] == 'FLUOROQUINOLONES'].loc[FQSL[FQSL['drug'] == 'FLUOROQUINOLONES']['percent of MDR isolates'].idxmax()]
	print("MAX measurable FQL is {}\t{}({}-{})".format(max_FQSL_FQ['country'], max_FQSL_FQ['percent of MDR isolates'], max_FQSL_FQ['lower_bound'], max_FQSL_FQ['upper_bound']))
	min_FQSL_FQ = FQSL[FQSL['drug'] == 'FLUOROQUINOLONES'].loc[FQSL[FQSL['drug'] == 'FLUOROQUINOLONES']['percent of MDR isolates'].idxmin()]
	print("MIN measurable FQL is {}\t{}({}-{})".format(min_FQSL_FQ['country'], min_FQSL_FQ['percent of MDR isolates'], min_FQSL_FQ['lower_bound'], min_FQSL_FQ['upper_bound']))
	max_FQSL_SL = FQSL[FQSL['drug'] == 'AMINOGLYCOSIDE'].loc[FQSL[FQSL['drug'] == 'AMINOGLYCOSIDE']['percent of MDR isolates'].idxmax()]
	print("MAX measurable SLI is {}\t{}({}-{})".format(max_FQSL_SL['country'], max_FQSL_SL['percent of MDR isolates'], max_FQSL_SL['lower_bound'], max_FQSL_SL['upper_bound']))
	min_FQSL_SL = FQSL[FQSL['drug'] == 'AMINOGLYCOSIDE'].loc[FQSL[FQSL['drug'] == 'AMINOGLYCOSIDE']['percent of MDR isolates'].idxmin()]
	print("MIN measurable SLI is {}\t{}({}-{})".format(min_FQSL_SL['country'], min_FQSL_SL['percent of MDR isolates'], min_FQSL_SL['lower_bound'], min_FQSL_SL['upper_bound']))


	#median MRSCA age for FLQ or SLI resistance acquisition among MDR isolates
	dates = dates[(dates['MDR'] == 1)].copy()
	# MFQSL = dates[(dates['FQ'] == 1) | (dates['SLIS'] == 1)]
	#No! the above is wrong we have to use the dates['drug'] since that is teh drug that was dated! other way could give us wrong dates
	MFQSL = dates[(dates['drug'] == 'FLUOROQUINOLONES') | (dates['drug'] == 'AMINOGLYCOSIDE')]
	print("MEDIAN MRSCA age for FLQ or SLI res aquisition among MDR isolates {}".format(MFQSL['date'].quantile([0.25,0.5,0.75])))


#Calculations for old_resistance stats in paragraph that starts with "The frequency of old resistance i.e. with MRSCA ≥20 years " 
def old_resistance():
	dates = pd.read_csv('strain_name_with_date.tsv',sep='\t')[['country','strain','drug','date']].copy()
	dates['drug'] = dates['drug'].replace({'PYRAZINAMIDE':'PZA','ISONIAZID':'INH','FLUOROQUINOLONES':'FLQ', 'STREPTOMYCIN':'STR',
		'RIFAMPICIN/RIFABUTIN':'RIF','AMINOGLYCOSIDE':'SLIs','ETHAMBUTOL':'EMB'})

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
	print("Percentage of old INH {} (95% CI {}-{})".format(oldINH*100, confidence_interval_lower(oldINHnum, oldINHdenom)*100, confidence_interval_upper(oldINHnum, oldINHdenom)*100))
	print("Percentage of old FLQ {} (95% CI {}-{})".format(oldFLQ*100, confidence_interval_lower(oldFLQnum, oldFLQdenom)*100, confidence_interval_upper(oldFLQnum, oldFLQdenom)*100))


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

old_resistance()

def geographic_distribution_MRSCA_comparison():
	dates = df.copy()
	dates = dates.replace({-1:0}).copy()
	to_keep = dates.groupby('country').sum().copy()
	del dates['date']
	name_with_date = pd.read_csv('strain_name_with_date.tsv',sep='\t')[['strain','drug','date']].copy()
	dates= pd.merge(dates, name_with_date, on='strain',how='outer')
	four_drug_interest = ['INH','RIF','SLIs','FLQ']
	five_countries_most_resistance = dates[(dates['S'] != 1)&(dates['country'] != 'Not Provided')&(dates['NO_DATA'] != 1)].groupby('country').count().nlargest(5, 'strain').reset_index()['country'].values
	print(five_countries_most_resistance)

	for drug in ['ISONIAZID', 'RIFAMPICIN/RIFABUTIN','AMINOGLYCOSIDE','FLUOROQUINOLONES','PYRAZINAMIDE','ETHAMBUTOL']:
		groups = []
		for country_one in five_countries_most_resistance:
			country_one_added = False
			for country_two in five_countries_most_resistance:
				if(country_one != country_two):
					one = dates[(dates['drug'] == drug)&(dates['country'] == country_one)]['date']
					two = dates[(dates['drug'] == drug)&(dates['country'] == country_two)]['date']
					# print("{}\t{}\t{}".format(country_one, country_two,sc.ranksums(one,two).pvalue))
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

	for drug in ['RIFAMPICIN/RIFABUTIN','AMINOGLYCOSIDE','FLUOROQUINOLONES']:
		print("IQR for United Kingdom {} is {}".format(drug, dates[(dates['drug'] == drug)&(dates['country'].isin(['United Kingdom']))]['date'].quantile([0.25,0.5,0.75])))

#Calculations for MRSCA ratio data in Supplementary Table 11
def MRSCA_ratio():
	name_with_date = pd.read_csv('strain_name_with_date.tsv',sep='\t')[['strain','drug','date']].copy()
	name_with_date['drug'] = name_with_date['drug'].replace({'PYRAZINAMIDE':'PZA','ISONIAZID':'INH','FLUOROQUINOLONES':'FLQ', 'STREPTOMYCIN':'STR',
		'RIFAMPICIN/RIFABUTIN':'RIF','AMINOGLYCOSIDE':'SLIs','ETHAMBUTOL':'EMB'})
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
	result = '('+test.round(3).groupby('country')['date_y'].nunique().astype('str')+'/'+country_total.astype('str')+')'+ (test.round(3).groupby('country')['date_y'].nunique()/country_total).round(2).astype('str')
	final['pooled'] = result.to_dict()


	####FOR SPECIFIC DRUGS#####

	for drug in ['INH','RIF','PZA','FLQ','STR','SLIs','EMB']:
		#How many dated strains in each country?
		country_total = test[test['drug'] == drug].groupby('country').sum()
		country_total = country_total[country_total['SUM'] > 10]
		country_total = country_total['SUM']
		#divide for answer
		result = '('+test[test['drug'] == drug].round(3).groupby('country')['date_y'].nunique().astype('str')+'/'+country_total.astype('str')+')'+ (test[test['drug'] == drug].round(3).groupby('country')['date_y'].nunique()/country_total).round(2).astype('str')
		final[drug] = result

	resultat = pd.DataFrame.from_dict(final).to_csv('MRSCA_ratios.tsv',sep='\t')
	print(resultat)



"""
Beneath this is old stuff that I keep in case we need it

"""


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






