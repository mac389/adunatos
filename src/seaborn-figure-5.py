import matplotlib, operator
matplotlib.use('Agg')
import string, json, brewer2mpl, mygene

import utils as tech
import seaborn as sns
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec

from awesome_print import ap 
from collections import Counter, OrderedDict

flattened_structural_ontology = open('../docs/flattened_structural_ontology','rb').read().splitlines()
df = pd.read_json('../docs/gene-expression-by-area2.json').dropna(axis=1)
genes_functions = json.load(open('../docs/gene-function.json','rb'))
cutoff = 3

#Direct input and output

#areas = tech.allChildrenOfParent('hippocampal formation',flattened_structural_ontology)
areas =open('tmp').read().splitlines()
ap(open('tmp').read().splitlines())
areas = [area for area in open('tmp').read().splitlines() if area in df.columns.values]
#ap(df.columns.values)
ap(areas)
data = df[areas]

#Create color palette for structures
regions = tech.array_from_lists([tech.ancestors(val,flattened_structural_ontology) for val in df.columns.values])
ind = np.lexsort(regions,axis=0)
regions = tech.nearest_not(regions[:,ind],-4,'None')
unique_regions = set(regions)

structure_color_palette = sns.color_palette('Set2',len(unique_regions))
structure_to_color_mapping = dict(zip(unique_regions,structure_color_palette))
col_colors = pd.Series(regions).map(structure_to_color_mapping)

#Create color palette for gene functions
functions = [genes_functions[gene] if type(genes_functions[gene])==str or type(genes_functions[gene]) == unicode else genes_functions[gene][-1] 
				for gene in df.index.values]
unique_functions = set(functions)

gene_color_palette = sns.diverging_palette(10, 220, sep=80, n=len(unique_functions))
gene_function_to_color_mapping = dict(zip(unique_functions,gene_color_palette))
row_colors = pd.Series(functions).map(gene_function_to_color_mapping)

#Determine frequency of each function
heatmapGS = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0,width_ratios=[0.5,1.5],height_ratios=[0.25,1])

function_frequences = OrderedDict(sorted(Counter(functions).items(),key=operator.itemgetter(1),reverse=True))
cg = sns.clustermap(data,method='complete', row_colors=row_colors,yticklabels=False,
	cbar_kws={'label':'Expression'}, clip_on=False, col_colors = col_colors, robust=True,
	linewidths=0)

legend_for = cg.ax_col_dendrogram.get_position()
for label in function_frequences.keys()[:10]:
    cg.ax_col_dendrogram.bar(0, 0, color=gene_function_to_color_mapping[label],
                            label=tech.format(label), linewidth=0)
cg.ax_col_dendrogram.legend(ncol=1,loc=(legend_for.x0 - .8, legend_for.y0))
cg.cax.set_position([.05, .5, .03, .1])

structure_frequencies = OrderedDict(sorted(Counter(regions).items(),key=operator.itemgetter(1),reverse=True))
for label in structure_frequencies.keys()[:10]:
    cg.ax_row_dendrogram.bar(0, 0, color=structure_to_color_mapping[label],
                            label=tech.format(label), linewidth=0)
legend_struct = cg.ax_row_dendrogram.get_position()
cg.ax_row_dendrogram.legend(ncol=1,loc=(legend_for.x0+2.5, legend_for.y0+0.65))
cg.savefig('SZ_hypo_test.png')
