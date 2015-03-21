import json, collections
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 

from copy import copy
from scipy.stats import zscore
from awesome_print import ap 
from jsonpath_rw import jsonpath, parse, Descendants
from scipy.stats.distributions import norm
from time import time
#plt.xkcd()

def before(x,y,full_string):
	return full_string.find(x) < full_string.find(y)

def isChildOfParent(child,parent,flattened_ontology):
	return len([path for path in flattened_ontology if before(parent,child,path)])>0

def traverse(list_of_dictionaries,_prefix=""):
	_items = []
	for structure in list_of_dictionaries:
		if structure['children'] != []:
			_items.extend(traverse(structure['children'],_prefix+'_'+structure['name']))
		else:
			_items.append(_prefix)
	return _items

def get_lineage(child,flattened_ontology):
	for lineage in flattened_ontology:
		print lineage

def common_ancestor(children,flattened_ontology):
	ap(children)
	lineages = [lineage.split('_') 
					for lineage in flattened_ontology 
					for child in children 
					if child in lineage]
	ap(lineages)
	'''				
	for child in children:
		for lineage in flattened_ontology:
			if child in lineage:
				ap('----------')
				ap(child)
				ap(lineage.split('_'))
				ap(lineage.split('_').index(child))
				ap('----------')
	'''
df = pd.read_json('./docs/gene-expression-by-area2.json').dropna(axis=1)
structural_ontolgy = json.load(open('./docs/brain-structure-ontology.json','rb'))['msg'][0]
#Group rows by region
#labels = {'Frontal Lobe':4009,'Insula':4268,'Limbic Lobe':4219,'Occipital Lobe':4180,
#			'Parietal Lobe':4084,'Temporal Lobe':4132} #White matter ventricles

labels = map(lambda s: '_%s_'%s,['frontal lobe','insula','limbic lobe','occipital lobe','parietal lobe','temporal lobe'])

flattened_structural_ontology = open('./docs/flattened_structural_ontology','rb').read().splitlines()

column_names = list(df.columns.values)
new_column_names = [name for name in column_names  for label in labels 
						if isChildOfParent(name,label,flattened_structural_ontology)]

new_column_names = list(collections.OrderedDict.fromkeys(new_column_names))

common_ancestor(new_column_names[:2],list(set(flattened_structural_ontology)))

df = df[new_column_names]
df = (df - df.mean())/df.std() #Expressions are in log value, so subtraction is scaling

'''
		Figure 4 structured as 

			Area
		|--------->
		|
Area    |    x_{ij}  represents the *number* of common genes enriched in one structure over another
		|                *still not sure what this phrase means*
		V

   ;Have to correct for simultaneous hypotheses; Bonferroni correction very conservative
   ; Didn't record variance across probes for a given gene, do know variance of gene expression across samples
   ; Using z-score assumes that the gene expression in each region follows a normal distribution (verify this)

'''

def compare(one,two):
	zscores = zscore(one-two)
	threshold = -norm.ppf(0.025/zscores.shape[0]) #I don't understand why there is a negative here
	return sum(np.absolute(zscores)>threshold)

differ = np.reshape([compare(df[col2],df[col]) for col2 in df for col in df], (df.shape[1],df.shape[1]))
differ = pd.DataFrame(differ, columns=df.columns.values,index=df.columns.values)

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(differ,cmap=plt.cm.jet,interpolation='nearest',aspect='auto')
cbar = plt.colorbar(cax)
cbar.set_label('No. of differentially expressed genes')
plt.savefig('./images/deg.tiff')