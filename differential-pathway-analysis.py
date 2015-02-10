import re

import pandas as pd 
import numpy as np 

from awesome_print import ap 

CUTOFF = 1

def parse_line(line):
	try:
		area,gene_level = re.split(r'[ ](?=[A-Z])',line)
		gene,level = gene_level.split()
		#area,gene,level = line.split('|')
		return {'Structure':area,'ID':gene,'Level':float(level)}
	except:
		return None


control = pd.DataFrame(filter(None,[parse_line(line) 
	for line in open('gene-expression-by-area-snapshot.txt').read().splitlines()]))

roi = pd.read_csv('gene-expression-by-area-filtered.txt',sep='|',
	names=['Structure','ID','Level'])

control['zscore'] = np.absolute((control['Level'] - control['Level'].mean())/control['Level'].std(ddof=0))
roi['zscore'] = np.absolute((roi['Level'] - roi['Level'].mean())/roi['Level'].std(ddof=0))


#Selectively upregulated genese are those in ROI with |zscore| >CUTOFF and with |zscore| < CUTOFF in Brain
roi_upregulated_genes = roi[roi['zscore']>CUTOFF]
control_not_upregulated_genes = control[control['zscore']<CUTOFF]

roi_selectively_upregulated_genes = pd.merge(roi_upregulated_genes,control_not_upregulated_genes,on='ID')
roi_selectively_upregulated_genes.drop_duplicates(cols='ID',take_last=True,inplace=True)

roi_selectively_upregulated_genes.to_pickle('differentially-expressed-genes.pkl')
