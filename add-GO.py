import json

import pandas as pd

from awesome_print import ap 
#This module adds the GO for each uniprot
csv = True
if csv:
	df = pd.read_csv('differentially-expressed-genes-with-uniprot.csv')
else:
	df = pd.read_pickle('differentially-expressed-genes-with-uniprot.pkl')

uniprot_to_goa = json.load(open('uniprot->goa_deduplicated.json','rb'))

go_list = [uniprot_to_goa[gene] if gene in uniprot_to_goa else ' '
			for gene in df['uniprot'].tolist()]

df['GO'] = go_list
df.to_pickle('differentially-expressed-genes-with-uniprot-go.pkl')
df.to_csv('differentially-expressed-genes-with-uniprot-go.csv')