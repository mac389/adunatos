import re, json, optparse,itertools
import utils as tech
from pprint import pprint
from awesome_print import ap 
from pprint import pprint 

broadman = json.load(open('broadman.json','rb'))
parser = optparse.OptionParser()
ontology = open('../docs/flattened_structural_ontology','rb').read().splitlines()

parser.add_option('--input','-i',
            action="store", dest="input",
            help="query string", default="")

options, args = parser.parse_args()

def parse_cluster(cluster,size_cutoff=1000):
	#Extract size
	data = re.split('(\d{2,5}\smm)',cluster)
	size = data[1]
	rest = data[2:]
	if size>size_cutoff:
		#Determine locations
		locs = rest[0].split('Labels:')[1:][0]
		locs = re.split('\d{2,5}mm\s',locs)[1:]
		areas = []
		for loc in locs:
			loc = loc.strip().split('.')[-1].split()[1:]
			key  = ' '.join(loc).capitalize()
			areas += [broadman[key]] if type(broadman[key]) == str or type(broadman[key]) == unicode else broadman[key] 
			#Entries are either strings or lists of strings
		areas = list(set(areas))
		return areas

text = open(options.input).read().split('Cluster Analysis:')[1:][0].split('Experiment Table')[:-1][0]
text = text.replace('\n','').strip()

areas = re.split('(#\d:)',text)[1:]
areas = [area.strip() for area in areas]
areas_dict = dict(zip(areas[::2],map(parse_cluster,areas[1::2])))

ontology_areas =list(itertools.chain.from_iterable(areas_dict.values()))
for area in ontology_areas:
	print area

