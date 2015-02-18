import re, itertools
import pandas as pd

from awesome_print import ap 

def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return itertools.izip(a, a)

da_genes_in_uniprot = pd.read_csv('../docs/human_da_signaling_genes',sep='\t')
da_genes_in_uniprot.columns = ['id_type','ID','ID2','pf?quoi','name','taxonomy','protein','GO','rien']

'''
data = re.split('(\[Term\]|\[Typdef\])',open('../docs/go-basic.obo','rb').read())

header = 1
footer = -1

parseable_data = list(pairwise(data[header:(footer-1)]))

big_go = {}
for term_type,content in parseable_data:
	if 'term' in term_type.lower():
		content_dict = {}
		for line in content.splitlines():
			if line != '':
				key,value = line.split(': ',1)
				if key not in content_dict:
					content_dict[key] = value
				else: #Can't use isinstance. Strings inherit from list.
					tmp = []
					tmp.append(content_dict[key])
					tmp.append(value)
					content_dict[key] = tmp

	big_go[content_dict['id']] = content_dict

ap(big_go['GO:0007212'])
ap(big_go["GO:0003674"])
ap([gene for id,gene in big_go.iteritems() if 'is_a' in gene and gene['is_a']=="GO:0007212"])
'''