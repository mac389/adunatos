import re 

import matplotlib.pyplot as plt 
import numpy as np 

from awesome_print import ap 
from matplotlib import rcParams
from scipy.signal import correlate2d

#plt.xkcd()
from matplotlib import rcParams
rcParams['text.usetex'] = True

cmap = plt.cm.seismic
cmap.set_bad('w',1.)
def parse_line(line):
	try:
		area,gene_level = re.split(r'[ ](?=[A-Z])',line)
		gene,level = gene_level.split()
		#area,gene,level = line.split('|')
		return {'area':area,'gene':gene,'level':float(level)}
	except:
		return None

def vcorrcoef(X):
	pass


data = filter(None,[parse_line(line) for line in open('gene-expression-by-area-snapshot.txt','rb').read().splitlines()])

areas,genes = zip(*[(item['area'],item['gene']) for item in data])
areas = list(set(areas))
genes = list(set(genes))

record = {}
for area in areas:
	record[area] = {}

for item in data:
	record[item['area']][item['gene']] = item['level']

expression_matrix = np.zeros((len(genes),len(areas)))

for row,gene in enumerate(genes):
	for col,area in enumerate(areas):
		if gene in record[area]:
			expression_matrix[row,col] = record[area][gene]
		else:
			expression_matrix[row,col] = np.nan


masked_expression_array = np.ma.array(expression_matrix,mask=np.isnan(expression_matrix))
#detrended_masked_expression_array = masked_expression_array - np.nanmean(masked_expression_array,axis=0)


ap('Creating expression matrix')

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(masked_expression_array, interpolation='nearest', cmap=cmap, aspect='auto')
ax.xaxis.tick_top()
ax.set_xticks([])
ax.set_yticks([])
#ax.set_xticklabels(areas,rotation='vertical')
ax.set_yticklabels([])
ax.set_ylabel(r'\Huge \textbf{\textsc{Genes}} $\rightarrow$')
ax.set_xlabel(r'\Huge \textbf{\textsc{Areas}} $\rightarrow$')
#cbar = plt.colorbar(cax)
#cbar.set_label('Expression Level')
plt.tight_layout()
plt.savefig('expression-overall-for-methods.png')
'''
del fig,ax 
ap('Creating correlation matrix')
'''

'''
  Corrcoeff expects       ------->
  						  |
  					Area  |
  						  |
  						  V

'''
'''
headers = np.loadtxt('gene-expression-by-area-filtered.txt',delimiter='|',dtype=str)
np.savetxt('gene-gene-correlation-subset.txt',corr_matrix,delimiter=',',
			header=' '.join(gene for gene in headers[:,1]))
'''
'''
corr_matrix = np.ma.corrcoef(masked_expression_array)
#corr_matrix = np.corrcoef(expression_matrix.T)
#mask = np.all(np.isnan(corr_matrix),axis=1)
#corr_matrix_cleansed = corr_matrix[~mask]
ap(corr_matrix.shape)
fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(corr_matrix, interpolation='nearest', cmap=cmap, aspect='equal',
	vmin=-1,vmax=1)
ax.xaxis.tick_top()
#ax.set_xticks(range(len(areas)))
#ax.set_xticklabels(areas,rotation='vertical')
#ax.set_yticks(range(len(areas)))
#ax.set_yticklabels(areas)
#ax.set_yticklabels([])
#ax.set_xticklabels([])
#ax.set_ylabel('Areas')
#ax.set_xlabel('Areas')
cbar = plt.colorbar(cax)
cbar.set_label('Correlation')
cbar.set_ticks([-1,0,1])
plt.tight_layout()
plt.savefig('area-area-correlation-for-methods-figure.png')

'''

