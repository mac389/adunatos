import mygene, json, os

import pandas as pd

'''
      This figure is structured as 
           Area  
      |----------->
Gene  |
      |   ij indicates the expression of gene i in area j
      |   height information is redundant. 
      V

'''

'''
 1. Get genes associated with DA signaling (later make this a command line argument)
 2. From GeneID -> Uniprot ID -> GO 
 2. Find expression by area
'''


gene_expression_by_area = pd.read_json('../gene-expression-by-area2.json').dropna(axis=1)
uniprot_to_goa = json.load(open('../docs/uniprot->goa_deduplicated.json','rb'))
mg = mygene.MyGeneInfo()
pathways_of_interest = {'DA':["GO:0001588","GO:0007212"],'GABA':['GO:0007214'],
  'Glu':['GO:0007215'],'Inflammation':["GO:0006954"],
  'Epigenetics':["GO:0040029"]}
gene_names = pd.DataFrame(list(set(list(gene_expression_by_area.index))), columns=['ID'])

if not os.path.isfile('gene_id_uniprogt_go_conversion.pkl'):
	kwargs = {'species':'human','fields':'uniprot','scopes':'symbol'}
	uniprot = []
	for i,gene in enumerate(gene_names['ID'].tolist()):
		txt = False
		if i%100==0:
			print 'On %d (%d)'%(i,len(gene_names['ID']))
		info = mg.query(gene,**kwargs)
		if 'hits' in info:
			info = info['hits'][0]
			if 'uniprot' in info and 'Swiss-Prot' in info['uniprot']:
				if type(info['uniprot']['Swiss-Prot']) == type([]):
					txt = info['uniprot']['Swiss-Prot'][0]
				else:
					txt = info['uniprot']['Swiss-Prot']
		uniprot.append(txt)
	gene_names['uniprot'] = uniprot

	go_list = [uniprot_to_goa[gene] if gene in uniprot_to_goa else ' '
				for gene in gene_names['uniprot'].tolist()]

	gene_names['GO'] = go_list

	gene_names.to_pickle('gene_id_uniprogt_go_conversion.pkl')
else:
	gene_names = pd.read_pickle('gene_id_uniprogt_go_conversion.pkl')

genes = gene_names['ID'][[i for i,go_id in enumerate(gene_names['GO']) 
				if len(set(pathways_of_interest['DA']) & set(go_id))>0]]


print len(genes)
	
'''
heatmap = np.array([[expression[gene][area] 
		for gene in genes] 
		for area in areas])
'''