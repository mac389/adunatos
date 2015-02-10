import mygene

import pandas as pd

from awesome_print import ap 
#This module identifies to which pathways the selectively upregulated genes belong
selectively_upregulated_genes = pd.read_pickle('differentially-expressed-genes.pkl')

mg = mygene.MyGeneInfo()


kwargs = {'species':'human','fields':'uniprot','scopes':'symbol'}
uniprot = []
for i,gene in enumerate(selectively_upregulated_genes['ID'].tolist()):
	txt = False
	if i%100==0:
		print 'On %d (%d)'%(i,len(selectively_upregulated_genes['ID']))
	info = mg.query(gene,**kwargs)
	if 'hits' in info:
		info = info['hits'][0]
		if 'uniprot' in info and 'Swiss-Prot' in info['uniprot']:
			if type(info['uniprot']['Swiss-Prot']) == type([]):
				txt = info['uniprot']['Swiss-Prot'][0]
			else:
				txt = info['uniprot']['Swiss-Prot']
	uniprot.append(txt)

selectively_upregulated_genes['uniprot'] = uniprot
selectively_upregulated_genes.to_pickle('differentially-expressed-genes-with-uniprot.pkl')
