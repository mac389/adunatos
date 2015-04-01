import re, json, optparse,itertools, os
import utils as tech
from pprint import pprint
from awesome_print import ap 
from pprint import pprint 

broadman = json.load(open('/Users/nirgal/Desktop/Frangou/Allen Brain Atlas/src/broadman.json','rb'))
parser = optparse.OptionParser()
ontology = open('/Users/nirgal/Desktop/Frangou/Allen Brain Atlas/docs/flattened_structural_ontology','rb').read().splitlines()

'''
parser.add_option('--input','-i',
            action="store", dest="input",
            help="query string", default="")

options, args = parser.parse_args()
'''

filenames = open('filelist').read().splitlines()

cutoff = 1000

for filename in filenames:
	outputname,_ = os.path.splitext(filename)
	condition = os.path.basename(outputname)
	outputname+='.structures'
	outputname = os.path.join(os.getcwd(),'ALE-output',outputname)
	with open(outputname,'wb') as outfile:
		text = open(os.path.join(os.getcwd(),'ALE-output',filename)).read().splitlines()
		#Find where clusters start
		text = ' '.join(text[(text.index('Cluster Analysis:')+1):text.index('Experiment Table:')])

		#split into clusters
		clusters = re.split('#\d{1,2}:',text)[1:]
		for cluster in clusters:
			#ap(cluster)
			cluster_size = re.split('(\d+\smm)',cluster)[1].split()[0]
			if int(cluster_size)> cutoff:
				structures = [structure.strip() for structure in re.findall(r'\b[^0-9()\*]+\s?\.[^0-9()\*]+\s?\w+',cluster)]
				structures = [structure.split('.') for structure in set(structures)]
				names = []
				for structure in structures:
					if any(['gyrus' in subfield.lower() for subfield in structure]):
						gyrus_idx = [structure.index(subfield) for subfield in structure if 'gyrus' in subfield.lower()].pop()
						names.append(structure[gyrus_idx])
					else:
						names.append(structure[-1])
				
				names = list(set([re.sub(r"(?<!area)\s{1}\d+","",re.sub(r"\d+mm"," ",name)).strip().lower() for name in names]))
				for idx,name in enumerate(names):
					if 'brodmann' in name:
						name = broadman[name.replace('brodmann','').strip().capitalize()]
					if type(name) == str or type(name) == unicode:
						print>>outfile,name
					if type(name) == list:
						for subfield in name:
							print>>outfile,subfield
	#Python equivalent of uniq from bash
	structures = set(open(outputname).read().splitlines())
	with open(outputname,'wb') as outfile:
		for structure in structures:
			print>>outfile,structure	

	ap('Structures extracted from %s'%condition)
