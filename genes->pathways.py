import json

from bioservices.kegg import KEGG
from awesome_print import ap 

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

ap(pathways)

'''
final_list = []
for item in pathways:
	if hasattr(item,'__iter__'):
		for key,value in item.iteritems():
			final_list.append('%s %s'%(key,value))
	else:
		final_list.append(item)

final_list = list(set(final_list))
ap(final_list)
'''