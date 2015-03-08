import mygene, json, os

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
from awesome_print import ap 

from mpl_toolkits.mplot3d.axes3d import Axes3D

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

gene_expression_by_area = pd.read_json('./docs/gene-expression-by-area2.json').dropna(axis=1)
uniprot_to_goa = json.load(open('./docs/uniprot->goa_deduplicated.json','rb'))
mg = mygene.MyGeneInfo()
pathways_of_interest = {'DA':["GO:0001588","GO:0007212"],'GABA':['GO:0007214'],
  'Glu':['GO:0007215'],'Inflammation':["GO:0006954"],
  'Epigenetics':["GO:0040029"]}

gene_names = pd.DataFrame(list(set(list(gene_expression_by_area.index))), columns=['ID'])

if not os.path.isfile('./docs/gene_id_uniprogt_go_conversion.pkl'):
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
	gene_names.to_pickle('./docs/gene_id_uniprogt_go_conversion.pkl')
else:
	gene_names = pd.read_pickle('./docs/gene_id_uniprogt_go_conversion.pkl')


da_genes_in_uniprot = pd.read_csv('./docs/human_da_signaling_genes',sep='\t')
da_genes_in_uniprot.columns = ['id_type','ID','ID2','pf?quoi','name','taxonomy','protein','GO','rien']

da_genes_in_uniprot_ids = list(set(da_genes_in_uniprot['ID'])) #retrieves only 31 unique genes, from 90 I guess otheres are isoforms?
areas = gene_expression_by_area.columns.values
idxs=list(gene_expression_by_area.index.values)

heatmap = np.array([[gene_expression_by_area[area][idxs.index(gene)] if gene in idxs else 0
		for gene in da_genes_in_uniprot_ids] 
		for area in areas])

y = np.array(range(heatmap.shape[0])).astype(int)
x = np.array(range(heatmap.shape[1])).astype(int)
X,Y = np.meshgrid(x,y)
z = np.array([heatmap[y,x] for y in Y for x in X]).astype(float)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cax = ax.plot_surface(X,Y,heatmap,rstride=1, cstride=1, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False)
cbar = plt.colorbar(cax)
cbar.set_label('Gene Expression Level')
ax.set_ylabel('Area')
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_xlabel('Genes')
plt.savefig('./Images/test-figure-2-3d.tiff')