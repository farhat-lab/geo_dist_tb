"""
Written by Yasha Ektefaie
June 2018

python3 geography.py <Strain identification table> <lineage table> <resistance data table>

Takes in the tables and outputs a dataframe in the form of a tsv file with every strain, lineage, and resistance data in one place
Then outputs a seperate dataframe in the form of a tsv file that displays country and each antibiotic with the number of strains in the country that 
showed resistivity or susceptbility to that antibiotic

"""
import sys
import csv
import pandas as pd
import re

arguments = sys.argv

#Take in the arguments and make it into a csv DictReader
strain_identification = csv.DictReader(open(arguments[1]), delimiter='\t')
lineages = csv.DictReader(open(arguments[2]), delimiter='\t')
resistance = csv.DictReader(open(arguments[3]), delimiter='\t')

geography = {}

strain_info = {}
strain_to_date = {}
run_to_xref = {}
queryrun_to_xref = {}

#print("STARTING OUR SEARCH")
antibiotics = resistance.fieldnames[1:]

#We first use the strain_identification to get a strain ID and the location associated with that ID
#NOTE: a lot of the locations were similar but needed to all be in one region EX: Saint Mary's Hospital and South Africa are both in South Africa
#So I had to put a lot of if statements to take care of the subtelties of the data

for row in strain_identification:
	key = None
	if(row['public_xref']):
		key = row['public_xref']
		#print("FOUND KEY: {}".format(key))
		strain_info[key] = []
	elif(row['internal_xref']):
		key = row['internal_xref']
		#print("FOUND KEY: {}".format(key))
		strain_info[key] = []

	runs = row['runs'].split(',')
	for run in runs:
		run_to_xref[run] = key

	query_run = row['query_run'].split(',')
	for run in query_run:
		queryrun_to_xref[run] = key

	place = None

	if('Peru' in key):
		geography['Peru'].append(key)
		place = 'Peru'
	elif('South Africa' in row['geo_loc_name'] or 'Kwa-Magwaza Hospital' in row['geo_loc_name'] or 'St' in row['geo_loc_name'] or 
		'Christ' in row['geo_loc_name'] or 'King' in row['geo_loc_name'] or 'M3' in row['geo_loc_name'] or 'Ethembeni' in row['geo_loc_name']
		or 'Durban' in row['geo_loc_name'] or 'Chwezi' in row['geo_loc_name']):
		if('South Africa' not in geography):
			geography['South Africa'] = [key]
		else:
			geography['South Africa'].append(key)
		place = 'South Africa'
	elif('India' in row['geo_loc_name']):
		if('India' not in geography):
			geography['India'] = [key]
		else:
			geography['India'].append(key)
		place = 'India'
	elif('Uganda' in row['geo_loc_name']):
		if('Uganda' not in geography):
			geography['Uganda'] = [key]
		else:
			geography['Uganda'].append(key)
		place = 'Uganda'
	elif('Guinea' in row['geo_loc_name']):
		if('Guinea' not in geography):
			geography['Guinea'] = [key]
		else:
			geography['Guinea'].append(key)
		place = 'Guinea'
	elif('Kenya' in row['geo_loc_name']):
		if('Kenya' not in geography):
			geography['Kenya'] = [key]
		else:
			geography['Kenya'].append(key)
		place = 'Kenya'
	elif('USA' in row['geo_loc_name'] or 'Goodwin' in row['geo_loc_name']):
		if('United States of America' not in geography):
			geography['United States of America'] = [key]
		else:
			geography['United States of America'].append(key)
		place = 'United States of America'
	elif('Scotland' in row['geo_loc_name'] or 'Dundee' in row['geo_loc_name']):
		if('Scotland' not in geography):
			geography['Scotland'] = [key]
		else:
			geography['Scotland'].append(key)
		place = 'Scotland'
	elif('Germany' in row['geo_loc_name']):
		if('Germany' not in geography):
			geography['Germany'] = [key]
		else:
			geography['Germany'].append(key)
		place = 'Germany'
	elif('Moldova' in row['geo_loc_name']):
		if('Moldova' not in geography):
			geography['Moldova'] = [key]
		else:
			geography['Moldova'].append(key)
		place = 'Moldova'
	elif('Thailand' in row['geo_loc_name']):
		if('Thailand' not in geography):
			geography['Thailand'] = [key]
		else:
			geography['Thailand'].append(key)
		place = 'Thailand'
	elif('Mali' in row['geo_loc_name']):
		if('Mali' not in geography):
			geography['Mali'] = [key]
		else:
			geography['Mali'].append(key)
		place = 'Mali'
	elif('Canada' in row['geo_loc_name'] or 'Catherine' in row['geo_loc_name']):
		if('Canada' not in geography):
			geography['Canada'] = [key]
		else:
			geography['Canada'].append(key)
		place = 'Canada'
	elif('Switzerland' in row['geo_loc_name'] or 'Siloah Clinic' in row['geo_loc_name']):
		if('Switzerland' not in geography):
			geography['Switzerland'] = [key]
		else:
			geography['Switzerland'].append(key)
		place = 'Switzerland'
	elif('UK' in row['geo_loc_name'] or 'United Kingdom' in row['geo_loc_name'] or 'Scotland' in row['geo_loc_name']  ):
		if('United Kingdom' not in geography):
			geography['United Kingdom'] = [key]
		else:
			geography['United Kingdom'].append(key)
		place = 'United Kingdom'
	elif('Tanzania' in row['geo_loc_name']):
		if('United Republic of Tanzania' not in geography):
			geography['United Republic of Tanzania'] = [key]
		else:
			geography['United Republic of Tanzania'].append(key)
		place = 'United Republic of Tanzania'
	elif(not row['geo_loc_name'] or 'missing' in row['geo_loc_name'] or 'not provided' in row['geo_loc_name']):
		if('Not Provided' not in geography):
			geography['Not Provided'] = [key]
		else:
			geography['Not Provided'].append(key)
		place = 'Not Provided'
	elif(row['geo_loc_name'] not in geography):
		#print("FOUND LOCATION: {}".format(row['geo_loc_name']))
		geography[row['geo_loc_name']] = [key]
		place = row['geo_loc_name']
	else:
		geography[row['geo_loc_name']].append(key)
		place = row['geo_loc_name']


	strain_info[key].append(place)

	#We get the date here and then we do some processing to make sure we just get just the year and nothing else from formats
	potential_date = row['collection_date']
	firstscreen = re.findall('[0-9][0-9][0-9][0-9]',potential_date)
	secondscreen = re.findall('-[0-9][0-9]', potential_date)
	if(firstscreen):
		strain_to_date[key] = firstscreen[0]
	elif(secondscreen):
		number = int(secondscreen[0][1:])
		if(number > 18):
			strain_to_date[key] = "19" + secondscreen[0][1:]
		else:
			strain_to_date[key] = "20" + secondscreen[0][1:]
	else:
		strain_to_date[key] = potential_date



#Next we break down the lineage of each strain--since the lineage has strain names that are the bioproject/biosample I check for those too
#Before I declare that the strain was not found

for lineage in lineages:
	correct_lineage = lineage['lineage']
	if("lineage1" in correct_lineage):
		correct_lineage = "lineage1"
	elif("lineage2" in correct_lineage):
		correct_lineage = "lineage2"
	elif("lineage3" in correct_lineage):
		correct_lineage = "lineage3"
	elif("lineage4" in correct_lineage):
		correct_lineage = "lineage4"
	elif("lineage5" in correct_lineage):
		correct_lineage = "lineage5"
	elif("lineage6" in correct_lineage):
		correct_lineage = "lineage6"
	elif("lineage7" in correct_lineage):
		correct_lineage = "lineage7"


	#For some reason lineage stuff uses ID and not xref -- so j put both made two try statements for both situations

	try:
		strain_info[lineage['ID']].append(correct_lineage)
	except KeyError:
		try:
			strain_info[run_to_xref[lineage['ID']]].append(correct_lineage)
		except KeyError:
			try:
				strain_info[queryrun_to_xref[lineage['ID']]].append(correct_lineage)
			except KeyError:
				#print("LINEAGE NO ID FOUND ERROR: {}".format(lineage['ID']))
				pass


	try:
		strain_info[lineage['xref']].append(correct_lineage)
	except KeyError:
		try:
			strain_info[run_to_xref[lineage['xref']]].append(correct_lineage)
		except KeyError:
			try:
				strain_info[queryrun_to_xref[lineage['xref']]].append(correct_lineage)
			except KeyError:
				#print("LINEAGE NO ID FOUND ERROR: {}".format(lineage['ID']))
				pass

#For the strains that don't have lineages I make sure to add an empty list in their place so that the analysis later on doesn't get messed up
for strain in strain_info:
	info = strain_info[strain]
	if(len(info) < 2):
		strain_info[strain].append([])

#Now we break down the resistance again checking if any keys are bioprojects/biosamples and then adding this to the strain info

for strain in resistance:
	fail = False
	currentkey = strain['Isolate'] 
	if(currentkey not in strain_info):
		if(currentkey in run_to_xref):
			currentkey = run_to_xref[currentkey]
		elif(currentkey in query_run):
			currentkey = query_run[currentkey]
		else:
			#print("RESISTANCE NO ID FOUND ERROR: {}".format(currentkey))
			fail=True

	resistant = []
	sensitive = []
	if(not fail):
		for antibiotic in antibiotics:
			if(strain[antibiotic] == 'R'):
				resistant.append(antibiotic)
			elif(strain[antibiotic] == 'S'):
				sensitive.append(antibiotic)

		strain_info[currentkey].append(resistant)
		strain_info[currentkey].append(sensitive)
		#print("STRAIN: {} RESISTANT: {} susceptable: {}".format(currentkey, resistant, sensitive))


#We really want to have a dataframe so that joining by country becomes really easy and in general dataframes are nice to work with
#So I create a list of dictionaries for each strain with the info in strain_info but in dictionary format for each strain
#To make sure this dataframe is uniform and not irregular I had to check for many edge cases which is where the if statements come in

list_of_strain_info_dictionaries = []

for strain in strain_info:
	sub_dict = {}
	info = strain_info[strain]
	sub_dict['strain'] = strain
	resfail = False
	susfail = False

	location = info[0]
	#print("STRAIN: {} INFO: {}".format(strain, info))
	sub_dict['country'] = location

	try:
		lineage = info[1]
		#print(lineage)
		if(lineage == "lineage1"):
			sub_dict['lineage1'] = 1
		else:
			sub_dict['lineage1'] = 0

		if(lineage == "lineage2"):
			sub_dict['lineage2'] = 1
		else:
			sub_dict['lineage2'] = 0

		if(lineage == "lineage3"):
			sub_dict['lineage3'] = 1
		else:
			sub_dict['lineage3'] = 0

		if(lineage == "lineage4"):
			sub_dict['lineage4'] = 1
		else:
			sub_dict['lineage4'] = 0

		if(lineage == "lineage5"):
			sub_dict['lineage5'] = 1
		else:
			sub_dict['lineage5'] = 0

		if(lineage == "lineage6"):
			sub_dict['lineage6'] = 1
		else:
			sub_dict['lineage6'] = 0

		if(lineage == "lineage7"):
			sub_dict['lineage7'] = 1
		else:
			sub_dict['lineage7'] = 0
	except IndexError:
		sub_dict['lineage1'] = 0
		sub_dict['lineage2'] = 0
		sub_dict['lineage3'] = 0
		sub_dict['lineage4'] = 0
		sub_dict['lineage5'] = 0
		sub_dict['lineage6'] = 0
		sub_dict['lineage7'] = 0

	try:
		resistant_antibiotics = info[2]
	except IndexError:
		resistant_antibiotics = []


	try:
		sensitive_antibiotics = info[3]
	except IndexError:
		sensitive_antibiotics = []

	#Logic if we have data for resistance it's resistant, if we have data for susceptible it's susceptible, else no data
	for resistant in antibiotics:
		if(resistant in resistant_antibiotics and resistant in sensitive_antibiotics):
			raise Exception("Something seriously went wrong this should not have happened: {}".format(strain))
		if(resistant in resistant_antibiotics):
			sub_dict[resistant] = 1
		elif(resistant in sensitive_antibiotics):
			sub_dict[resistant] = 0
		else:
			sub_dict[resistant] = -1

	# if(not susfail):
	# 	for susceptable in antibiotics:
	# 		if(susceptable in sensitive_antibiotics):
	# 			sub_dict[susceptable+"_S"] = 1
	# 		else:
	# 			sub_dict[susceptable+"_S"] = 0

	#Just really quickly adding in the date here
	sub_dict['date'] = strain_to_date[strain]

	list_of_strain_info_dictionaries.append(sub_dict)

#With the list of dictionaries I can now make the dataframe of interest, copy it, delete the strain and lineage parts, and then join on the country and sum!
strain_info_df = pd.DataFrame(list_of_strain_info_dictionaries)

#A quick aside where I create a dataframe mapping country to number of strains in that country--this is needed for some data analysis
list_of_geography = []
for key in geography:
	new_geography = {}
	new_geography['country'] = key
	new_geography['number'] = len(geography[key])
	list_of_geography.append(new_geography)

country_numbers = pd.DataFrame(list_of_geography)

#NEW STUFF ADDED IN####
strains = [i.rstrip() for i in open('final_all_strains','r').readlines()]
strain_info_df = strain_info_df[strain_info_df['strain'].isin(strains)]


#I save the country numbers mapping info and the strain info dictionary
strain_info_df.to_csv('strain_info.tsv', sep='\t')
country_numbers.to_csv('country_numbers.tsv',sep='\t')

#I do the copying, deletion, and summing I talked about before here
geo_resistance_data = strain_info_df.copy(deep=True)
del geo_resistance_data['strain']

country_resistance = geo_resistance_data.groupby("country").sum()

#Now I have to normalize all of my data against the number of strains in each country--I do that here
indices_of_interest = [i+"_R" for i in antibiotics]+[i+"_S" for i in antibiotics]+['lineage1','lineage2', 'lineage3', 'lineage4']

for i, trial in country_resistance.iterrows():
	number = len(geography[trial.name])
	#print("COUNTRY: {} NUMBER: {}".format(trial.name, number))
	country_resistance.loc[i] = country_resistance.loc[i]/number

#I save this normalized, cut, modified data here
country_resistance.to_csv('resistance_info.tsv', sep='\t')
