import itertools
import brewer2mpl
import sys

import Graphics as artist
import pandas as pd 
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.spatial.distance import pdist, squareform
from scipy.cluster import hierarchy 
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib.colors import LinearSegmentedColormap
from awesome_print import ap 
from matplotlib import rcParams
from optparse import OptionParser

sys.setrecursionlimit(10000)

parser = OptionParser(usage="usage: %prog [options] filename",
													version="%prog 1.0")
parser.add_option("-s", "--structures",
									action="store",
									dest="target_area",
									help="structure name") #Later expand this to also take a file name 

parser.add_option("-o", "--output",
									action="store",
									dest="savename",
									help="where to save images")

options, args = parser.parse_args()
print options.target_area
# IN FIGURE 5 THEY USE DIFFERENTIALLY EXPRESSED GENES, NOT EXPRESSION LEVEL, THIS FIGURE IS USING EXPRESSION LEVEL, CHANGE IT. 

#Abstract by making the area something you input

rcParams['font.size'] = 10
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

cdict = {'red':   ((0.0, 0.0, 0.0),
									 (0.5, 0.0, 0.0),
									 (1.0, 1.0, 1.0)),
				 'blue':  ((0.0, 0.0, 0.0),
									 (1.0, 0.0, 0.0)),
				 'green': ((0.0, 0.0, 1.0),
									 (0.5, 0.0, 0.0),
									 (1.0, 0.0, 0.0))}

cmap = LinearSegmentedColormap('nature', cdict, 100)

def clean_axis(ax):
		"""Remove ticks, tick labels, and frame from axis"""
		ax.get_xaxis().set_ticks([])
		ax.get_yaxis().set_ticks([])
		for sp in ax.spines.values():
				sp.set_visible(False)

hierarchy.set_link_color_palette(['black'])
structure_colorbar = brewer2mpl.get_map('Set2','qualitative',7)
gene_colorbar = brewer2mpl.get_map('Set3','qualitative',10)
lcolor = len(structure_colorbar.mpl_colors)
lgcolor = len(gene_colorbar.mpl_colors)
norm = mpl.colors.Normalize(vmin=5, vmax=10)
flattened_structural_ontology = open('./docs/flattened_structural_ontology','rb').read().splitlines()
df = pd.read_json('./docs/gene-expression-by-area2.json').dropna(axis=1)

'''
		 Data formatted as 
					Area
			|----------->
			|
Genes |   
			|
			V

'''

def after(path,parent):
	try:
		return [structure for structure in path if path.index(structure) > path.index(parent)]
	except:
		return [parent] #terminal structure

def before(path,parent):
	return [structure for structure in path if path.index(structure) < path.index(parent)]

def allChildrenOfParent(parent,ontology):
	paths = list(set(itertools.chain.from_iterable([after(path.split('_'),parent) 
				for path in ontology if "%s_"%parent in path or "_%s_"%parent in path])))
	if paths == []:
		return [parent] #terminal structure
	else:
		return paths

def ancestor(currentNode,ontology,levels=1): #Count of levels of ancestors to return
	paths = [path for path in ontology if "%s_"%currentNode in path or "_s_"%currentNode in path]
	ans = list(set(itertools.chain.from_iterable([[structure for structure in path.split('_')  
				if path.split('_').index(currentNode) > path.split('_').index(structure)] for path in paths])))[-levels:]
	#Flattening won't work structure paths of very different lengths, I don't think . 
areas = allChildrenOfParent(options.target_area,flattened_structural_ontology)

ap(areas)
areas = [area for area in areas if area in df.columns.values]

data = df[areas]#.tail(50)

col_labels = data.columns.values
row_labels = data.index.values

row_dist = pd.DataFrame(squareform(pdist(data,metric='euclidean')),
	columns=row_labels,index=row_labels) 

row_clusters = linkage(row_dist,method='complete')

cluster_df = pd.DataFrame(row_clusters,
		columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
	index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])

# Compute pairwise distances for columns
col_dists = squareform(pdist(df.T, metric='euclidean'))
col_clusters = linkage(col_dists, method='complete')
 # makes dendrogram black)
# reorder columns and rows with respect to the clustering

##$%^&*()
fig = plt.figure(figsize=(12,8))
heatmapGS = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0,width_ratios=[0.25,1],height_ratios=[0.25,1])

### row dendrogram ###
rowGSSS = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=heatmapGS[1,0],wspace=0.0,hspace=0.0,width_ratios=[1,0.25])
row_denAX = fig.add_subplot(rowGSSS[0,0])
row_dendrogram = dendrogram(row_clusters, orientation='right',color_threshold=np.inf) # makes dendrogram black
#ap(row_dendrogram)
df_rowclust = df.ix[row_dendrogram['leaves']]
clean_axis(row_denAX)

row_cbAX = fig.add_subplot(rowGSSS[0,1])
threshold=.2
row_cbSE = pd.Series([gene_colorbar.mpl_colors[i%lgcolor] for i in fcluster(row_clusters,threshold)])
row_axi = row_cbAX.imshow([ [x] for x in row_cbSE.ix[row_dendrogram['leaves']].values ],interpolation='nearest',aspect='auto',origin='lower')
clean_axis(row_cbAX)
row_colorbar_cbar = mpl.colorbar.ColorbarBase(fig.add_axes([0.04, 0.2, 0.03, 0.15]), 
							cmap=plt.cm.Set3, orientation='vertical', norm=norm)
row_colorbar_cbar.set_label('Genes', labelpad=-50)
row_colorbar_cbar.ax.set_yticklabels(['A','B','C','D','E',"F","G"])
row_colorbar_cbar.outline.set_visible(False)

### col dendrogram ####
colGSSS = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=heatmapGS[0,1],wspace=0.0,hspace=0.0,height_ratios=[1,0.25])
col_denAX = fig.add_subplot(colGSSS[0,0])
col_dendr = dendrogram(col_clusters,color_threshold=np.inf)
df_colrowclust = df_rowclust.ix[:][col_dendr['leaves']]
clean_axis(col_denAX)

### col colorbar ###
col_cbAX = fig.add_subplot(colGSSS[1,0])
threshold= 1
col_cbSE = pd.Series([structure_colorbar.mpl_colors[i%lcolor] for i in fcluster(col_clusters,threshold)])
col_axi = col_cbAX.imshow([list(col_cbSE.ix[col_dendr['leaves']])],interpolation='nearest',aspect='auto',origin='lower')
clean_axis(col_cbAX)
col_colorbar_cbar = mpl.colorbar.ColorbarBase(fig.add_axes([0.04, 0.02, 0.03, 0.15]), 
							cmap=plt.cm.Set1, orientation='vertical', norm=norm)
col_colorbar_cbar.set_label('Structures',labelpad=-50)
col_colorbar_cbar.ax.set_yticklabels(['A','B','C','D','E',"F","G"])
col_colorbar_cbar.outline.set_visible(False)

### heatmap ###
heatmapAX = fig.add_subplot(heatmapGS[1,1])
axi = heatmapAX.imshow(df_colrowclust,interpolation='nearest',aspect='auto',origin='lower',cmap=cmap)
clean_axis(heatmapAX)

## row labels ##
heatmapAX.set_yticks(np.arange(df_colrowclust.shape[0]))
heatmapAX.yaxis.set_ticks_position('right')
heatmapAX.set_yticklabels(df_colrowclust.index.values, fontsize=8)

## col labels ##
heatmapAX.set_xticks(np.arange(df_colrowclust.shape[1]))
xlabelsL = heatmapAX.set_xticklabels(df_colrowclust.columns.values, fontsize=10)
# rotate labels 90 degrees
for label in xlabelsL:
		label.set_rotation(90)
# remove the tick lines
for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines(): 
		l.set_markersize(0)

### scale colorbar ###
scale_cbGSSS = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=heatmapGS[0,0],wspace=0.0,hspace=0.0)
scale_cbAX = fig.add_subplot(scale_cbGSSS[0,1]) # colorbar for scale in upper left corner
cb = fig.colorbar(axi,scale_cbAX) # note that we could pass the norm explicitly with norm=my_norm
cb.set_label('Gene Expression Level')
cb.ax.yaxis.set_ticks_position('left') # move ticks to left side of colorbar to avoid problems with tight_layout
cb.ax.yaxis.set_label_position('left') # move label to left side of colorbar to avoid problems with tight_layout
cb.outline.set_linewidth(0)

# make colorbar labels smaller
tickL = cb.ax.yaxis.get_ticklabels()
for t in tickL:
		t.set_fontsize(t.get_fontsize() - 3)

heatmapGS.tight_layout(fig,h_pad=0.1,w_pad=0.5)
plt.savefig('./Images/%s.tiff'%options.savename 
		if not options.savename.endswith('.tiff') else options.savename)
