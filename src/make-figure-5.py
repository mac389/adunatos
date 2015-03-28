import itertools, brewer2mpl, mygene

'''
      CLEAN THIS CODE UP. IT WILL BE INCOMPREHENSIBLE TO OTHERS AND ME AFTER A HIATUS
'''

import Graphics as artist
import pandas as pd 
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import utils as tech

from scipy.spatial.distance import pdist, squareform
from scipy.cluster import hierarchy 
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib.colors import LinearSegmentedColormap
from collections import OrderedDict
from awesome_print import ap 
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator, FixedFormatter
# IN FIGURE 5 THEY USE DIFFERENTIALLY EXPRESSED GENES, NOT EXPRESSION LEVEL, THIS FIGURE IS USING EXPRESSION LEVEL, CHANGE IT. 

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

hierarchy.set_link_color_palette(['black'])
structure_colorbar = brewer2mpl.get_map('Set2','qualitative',7)
gene_colorbar = brewer2mpl.get_map('Set3','qualitative',10)
lcolor = len(structure_colorbar.mpl_colors)
lgcolor = len(gene_colorbar.mpl_colors)
norm = mpl.colors.Normalize(vmin=5, vmax=10)
flattened_structural_ontology = open('../docs/flattened_structural_ontology','rb').read().splitlines()
df = pd.read_json('../docs/gene-expression-by-area2.json').dropna(axis=1)
cutoff = 3

'''
     Data formatted as 
          Area
      |----------->
      |
Genes |
      |
      V

'''

areas = tech.allChildrenOfParent('hippocampal formation',flattened_structural_ontology)
areas = [area for area in areas if area in df.columns.values]

data = df[areas].tail(100)

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

fig = plt.figure(figsize=(20,8))
heatmapGS = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0,width_ratios=[0.5,1.5],height_ratios=[0.25,1])

### row dendrogram ###
rowGSSS = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=heatmapGS[1,0],wspace=0.0,hspace=0.0,width_ratios=[1,0.25])
row_denAX = fig.add_subplot(rowGSSS[0,0])
row_dendrogram = dendrogram(row_clusters, orientation='right',color_threshold=np.inf,ax=row_denAX) 
df_rowclust = df.ix[row_dendrogram['leaves']]
### col dendrogram ####
colGSSS = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=heatmapGS[0,1],
                                  wspace=0.0,hspace=0.0,height_ratios=[1,0.25])
col_denAX = fig.add_subplot(colGSSS[0,0])
col_dendr = dendrogram(col_clusters,color_threshold=np.inf,ax=col_denAX)
df_colrowclust = df_rowclust.ix[:][col_dendr['leaves']]
tech.clean_axis(row_denAX)
row_cbAX = fig.add_subplot(rowGSSS[0,1])
threshold=.2

ylabels = map(tech.get_gene_function,df_colrowclust.index.values)
idxs = tech.argsort(ylabels)
ylabels = [ylabels[i] for i in idxs]
gene_functions = {gene_function:gene_colorbar.mpl_colors[i%lcolor] 
                for i,gene_function in enumerate(ylabels)}
lineticks = tech.get_concept_boundaries(ylabels,cutoff=cutoff)
ylabels,yticks = zip(*tech.word_change_boundaries(ylabels,cutoff=cutoff).items())

gene_functions = OrderedDict(gene_functions.items(),key=lambda item:item[0])
row_cbSE = pd.Series([gene_functions[tech.get_proper_gene_function_key(df_colrowclust.index.values[i],gene_functions.keys())] 
                  for i,_ in enumerate(fcluster(row_clusters,threshold))])

row_axi = row_cbAX.imshow([[x] for x in row_cbSE.ix[row_dendrogram['leaves']].values ],
                                    interpolation='nearest',aspect='auto',origin='lower')
tech.clean_axis(row_cbAX)
row_colorbar_cbar = mpl.colorbar.ColorbarBase(fig.add_axes([0.16, 0.01, 0.03, 0.15]), 
							cmap=plt.cm.Set3, orientation='vertical', norm=norm)
row_colorbar_cbar.set_label('Functions', labelpad=-100)
row_colorbar_cbar.outline.set_visible(False)

tech.clean_axis(col_denAX)

labels = tech.array_from_lists([tech.ancestors(val,flattened_structural_ontology) 
                for val in df_colrowclust.columns.values])
ind = np.lexsort(labels,axis=0)
xlabels = tech.nearest_not(labels[:,ind],-4,'None')
#idxs = tech.argsort(xlabels)
#xlabels=[xlabels[idx] for idx in idxs]
xlineticks = tech.get_concept_boundaries(xlabels,cutoff=cutoff)
xlabels,xticklocs = zip(*tech.word_change_boundaries(xlabels,cutoff=cutoff).items())

#map structure to color
structures = {structure:structure_colorbar.mpl_colors[i%lcolor] 
                  for i,structure in enumerate(xlabels)}

structures = OrderedDict(structures.items(),key=lambda item: item[0])

### col colorbar ###
col_cbAX = fig.add_subplot(colGSSS[1,0])
threshold= 1

col_cbSE = pd.Series([structures[tech.get_proper_key(df_colrowclust.columns.values[i],structures,
            flattened_structural_ontology, fill_value='frontal lobe')] 
                for i,_ in enumerate(fcluster(col_clusters,threshold))])
#Do I need to enumerate over the clusters, or just the indices?
col_axi = col_cbAX.imshow([list(col_cbSE.ix[col_dendr['leaves']])],
        interpolation='nearest',aspect='auto',origin='lower')
tech.clean_axis(col_cbAX)
col_colorbar_cbar = mpl.colorbar.ColorbarBase(fig.add_axes([0.04, 0.01, 0.03, 0.15]), 
							cmap=plt.cm.Set1, orientation='vertical', norm=norm)
col_colorbar_cbar.set_label('Structures',labelpad=-150)
col_colorbar_cbar.ax.set_yticklabels([item[0].capitalize() for item in structures.items()]) 
col_colorbar_cbar.outline.set_visible(False)

### heatmap ###
heatmapAX = fig.add_subplot(heatmapGS[1,1])
axi = heatmapAX.imshow(df_colrowclust,interpolation='nearest',aspect='auto',origin='lower',cmap=cmap)
tech.clean_axis(heatmapAX)

heatmapAX.set_yticks(np.arange(df_colrowclust.shape[0]))
heatmapAX.yaxis.set_ticks_position('right')
ylabels = map(tech.get_gene_function,df_colrowclust.index.values)
idxs = tech.argsort(ylabels)
ylabels = [ylabels[i] for i in idxs]
lineticks = tech.get_concept_boundaries(ylabels,cutoff=cutoff)
ylabels,yticks = zip(*tech.word_change_boundaries(ylabels,cutoff=cutoff).items())
heatmapAX.set_yticks(yticks)
heatmapAX.set_yticklabels(ylabels, fontsize=8)
for linetick in lineticks:
  heatmapAX.axhline(linetick-1,color='k',clip_on=False, xmax=len(df_colrowclust.columns))
row_colorbar_cbar.ax.set_yticklabels([item[0].capitalize() for item in gene_functions.items()])

## col labels ##
heatmapAX.set_xticks(np.arange(df_colrowclust.shape[1]))
xlabelsL = heatmapAX.set_xticklabels(xlabels, fontsize=10)
for linetick in xlineticks:
  heatmapAX.axvline(linetick-1,color='k',clip_on=False, ymin=-len(df_colrowclust.columns))
  
## row labels ##
# rotate labels 90 degrees
for label in xlabelsL:
    label.set_rotation(90)
# remove the tick lines
for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines(): 
    l.set_markersize(0)
heatmapAX.xaxis.set_major_locator(FixedLocator(xticklocs))
heatmapAX.xaxis.set_major_formatter(FixedFormatter(xlabels))
### scale colorbar ###
'''
  Alter this code. Colorbar is stretched out horizontally. 
'''
scale_cbGSSS = gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=heatmapGS[0,0],wspace=0.0,hspace=0.0)
scale_cbAX = fig.add_subplot(scale_cbGSSS[0,1]) # colorbar for scale in upper left corner
cb = fig.colorbar(axi,scale_cbAX) # ould pass the norm explicitly with norm=my_norm
cb.set_label('Expression')
cb.ax.yaxis.set_ticks_position('left') # move ticks to left side of colorbar to avoid problems with tight_layout
cb.ax.yaxis.set_label_position('left') # move label to left side of colorbar to avoid problems with tight_layout
cb.outline.set_linewidth(0)

# make colorbar labels smaller
tickL = cb.ax.yaxis.get_ticklabels()
for t in tickL:
    t.set_fontsize(t.get_fontsize() - 3)

heatmapGS.tight_layout(fig,h_pad=0.1,w_pad=0.5)
plt.savefig('fig-5-small.tiff')
