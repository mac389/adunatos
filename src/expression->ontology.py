import json, mygene

import numpy as np
import pandas as pd

import goatools

from awesome_print import ap 
from Bio import Entrez
from bioservices.kegg import KEGG

Entrez.email = 'mac389@gmail.com'

association = json.load(open('go_association.json','rb'))
gene_ids_to_names = json.load(open('gene_ids','rb'))
uniprot_ids_to_go = json.load(open('uniprot->goa_deduplicated.json','rb'))
mg = mygene.MyGeneInfo()

def retrieve_annotation(id_list):
	request = Entrez.epost('gene',id=','.join(id_list))
	try:
		result = Entrez.read(request)
	except:
		print 'Error'

	webEnv = result["WebEnv"]
	queryKey = result["QueryKey"]
	data = Entrez.esummary(db="gene", webenv=webEnv, 
		query_key = queryKey)
	annotations = Entrez.read(data)
	return annotations


df = pd.read_csv('gene-expression-by-area-filtered.txt',sep='|',
	names=['Structure','ID','Level'])

by_level = df.sort(['Level'],ascending=False)
ids = list(by_level[:100]['ID'].values)

with open('gwas_gids.txt','wb') as f:
	for id in ids:
		if id in gene_ids_to_names:
			print>>f,gene_ids_to_names[id]

gene_ids = list(set(open('gwas_gids.txt','rb').read().splitlines()))
gene_ids_to_name = json.load(open('gene_ids_to_names.json','rb'))
k = KEGG()

pathways = {gene_ids_to_name[gene_id]:k.get_pathway_by_gene(gene_id,"hsa")
			for gene_id in gene_ids}

for gene in pathways:
	if type(pathways[gene]) == type({}):
		tmp = pathways[gene]
		tmp = {key:value for key,value in pathways[gene].iteritems() 
				if 'hsa' in key}
		pathways[gene] = tmp

json.dump(pathways,open('enriched-ontology.json','wb'))

with open('enriched-ontology.txt','wb') as outfile:
	for gene,pathway in pathways.iteritems():
		print>>outfile,'-----------'
		print>>outfile,gene
		print>>outfile,pathway
		print>>outfile,'-----------'



#To get GOA, convert Entrez ID to Uniprot ids

'''
uniprots = [gene['uniprot']['Swiss-Prot'] if type(gene['uniprot']['Swiss-Prot']) != type([])
		else gene['uniprot']['Swiss-Prot'][0] 
		for gene in mg.querymany(ids,scopes='symbol',fields='uniprot',species='human')]


overexpressed_ontologies = [uniprot_ids_to_go[uniprot_id] for uniprot_id in 
									uniprots]

ap(overexpressed_ontologies)
'''