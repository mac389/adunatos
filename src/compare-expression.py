import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import Graphics as artist

from awesome_print import ap 
from matplotlib import rcParams

rcParams['text.usetex'] = True
#This module compares expression between desired regions and the overall brain

def gene_in_pathway(gene_df,pathway_go):
	return len(set(gene_df['GO']) & set(pathway_go))>0

pathways_of_interest = {'DA':["GO:0001588","GO:0007212"],'GABA':['GO:0007214'],
  'Glu':['GO:0007215'],'Inflammation':["GO:0006954"],
  'Epigenetics':["GO:0040029"]}

df = pd.read_pickle('../docs/differentially-expressed-genes-with-uniprot-go.pkl')
sideviews = {}
for_frangou = {}
gene_gos = map(set,df['GO'].values)
for pathway,go_ids in pathways_of_interest.iteritems():
	idx = [len(gene_go & set(go_ids))>0 for gene_go in gene_gos]
	sideviews[pathway] = df[idx]
	#create views onto df 
	df[idx][['ID','Level_x','Level_y']].to_csv('../docs/differentially-expressed-genes-for-%s.txt'%pathway,
		index=False,header=["Gene", "Expression in ROI","Expression in Brain"])

#Visualization
ind = np.arange(len(pathways_of_interest))
pathway_keys = pathways_of_interest.keys()
width=0.35
fig = plt.figure()
ax = fig.add_subplot(111)
brain_mu = [sideviews[pathway]['zscore_y'].mean() for pathway in pathway_keys]
brain_std = [sideviews[pathway]['zscore_y'].std() for pathway in pathway_keys]

roi_mu = [sideviews[pathway]['zscore_x'].mean() for pathway in pathway_keys]
roi_std = [sideviews[pathway]['zscore_x'].std() for pathway in pathway_keys]

brain_patch = ax.bar(ind,brain_mu,width,color='k',alpha=0.8,yerr=brain_std,
	label=artist.format('Brain'),error_kw={'ecolor':'k','lw':2},edgecolor='k',clip_on=False)
roi_path = ax.bar(ind+width,roi_mu,width,color='r',alpha=0.8,yerr=roi_std,
	label=artist.format('ROI'),error_kw={'ecolor':'r','lw':2},edgecolor='r',clip_on=False)

artist.adjust_spines(ax)
ax.set_ylabel(artist.format('Z-score, Gene Expression Level'))
ax.set_xticks(ind+width)
ax.set_xticklabels(map(artist.format,pathway_keys))
for i,pathway in enumerate(pathway_keys):
	if i==2 or i==4:
		i+=.1
	ax.annotate(r'\textbf{\textsc{(%d)}}'%(len(sideviews[pathway]['zscore_x'])), xycoords='figure fraction',
		xy=(0.17+0.15*i,0.01))
plt.legend(frameon=False)

#json.dump(to_plot,open('for-frangou-presentation.json','wb'))

plt.savefig('../Images/enriched-genes-by-area-2.tiff')