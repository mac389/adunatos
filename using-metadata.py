import csv, requests, json, mygene

import numpy as np 

from awesome_print import ap 


def make_payload(structure):
	#return {"criteria":"model::Structure,rma::criteria[name$il'%s'],ontology[name$eq'Human Brain Atlas'],rma::options[only$eq'structures.id']"%structure}
	#return {"criteria":"model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq'%s'],rma::options[only$eq'probes.id']"%structure}
	return {"criteria":"model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[id$eq'%s'],rma::options[only$eq'probes.id']"%structure}

def make_full_payload(probe,donor,structure):
	#probes is a common-separated string of probe ids with no spaces
	return {"criteria":"service::human_microarray_expression[probes$eq%s][donors$eq%s][structures$eq%s]"%(probe,donor,structure)}

donors = list(set([item['donor_id'] for item in csv.DictReader(open('Columns.csv','rb'))]))
structures = json.load(open('structure_ids.json','rb'))
genes_to_probes = json.load(open('genes->probes.json','rb'))
gene_ids_to_names = json.load(open('gene_ids_to_names.json','rb'))
base ='http://api.brain-map.org'
query_url = base + '/' + 'api/v2/data/query.json'

donor = ','.join(map(str,donors))
genes_by_area = {}

total_iters = len(structures)*len(genes_to_probes)
iter_count = 0
tries = 0
excepts = 0
with open('gene-expression-by-area.txt','wb') as outfile:
	for structure_name,structure_id in structures.iteritems():
		genes_by_area[structure_name] = {}
		for gene_id,probes in genes_to_probes.iteritems():
			try:
				full_request = requests.get(query_url,params=make_full_payload(','.join(map(str,probes)),donor,structure_id))
				average_expression_level = np.array([response['expression_level'] 
					for response in full_request.json()['msg']['probes']]).astype(float).mean()
			
				genes_by_area[structure_name][gene_ids_to_names[str(gene_id)]] = average_expression_level
				print>>outfile,'%s %s %.04f'%(structure_name,gene_ids_to_names[str(gene_id)],average_expression_level)
				tries +=1
			except:
				excepts += 1
			iter_count += 1
			if iter_count%100==0:
				print '%d / %d (%.02f); tries =  %d | excepts = %d'%(iter_count,total_iters, iter_count/float(total_iters), tries,excepts)

json.dump(genes_by_area,open('gene-expression-by-area2.json','wb'))
'''
					    Expression
					------------ > 
					|
	Structure		|
					|
					V

			\
			/
			\
			|
			V
	PCA of the correlation matrix of above 
	MUST THINK ABOUT SCALING VARIANCES

'''