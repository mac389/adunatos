import pandas as pd
import matplotlib.pyplot as plt 

from scipy.stats import ranksums
from awesome_print import ap 

plt.xkcd()
df = pd.read_json('../gene-expression-by-area2.json').dropna(axis=1)

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
	column_one_header,column_one_data = one
	column_two_header,column_two_data = two

	#ap('%s -vs- %s'%(column_one_header,column_two_header))

	#Ranksums gives z and then p-value
	return ranksums(column_one_data,column_two_data)[0]

heatmap = [[compare(column_one,column_two)
			for column_one in df.iteritems()]
			for column_two in df.iteritems()]

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.imshow(heatmap,cmap=plt.cm.seismic,interpolation='nearest',aspect='auto')
cbar = plt.colorbar(cax)
cbar.set_label('Z-score; Differential Gene Expression')
plt.savefig('../images/for-nature-differential-gene-expression.png')
