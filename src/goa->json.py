import json
from awesome_print import ap 
from collections import defaultdict 

def process_line(line):
	tokens = line.split()
	return {tokens[1]:[token for token in tokens if 'GO:' in token and '_REF' not in token]}
data = open('gene_association.goa_human','rb').read().splitlines()

data = map(process_line,data)
#json.dump(map(process_line,data),open('uniprot->goa.json','wb'))

better_data = defaultdict(list)
for d in data:
	for key,value in d.iteritems():
		better_data[key].extend(value)

json.dump(better_data,open('uniprot->goa_deduplicated.json','wb'))