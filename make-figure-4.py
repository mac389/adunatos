import pandas as pd

df = pd.read_json('gene-expression-by-area2.json').dropna(axis=1)

'''
        Figure 4 structured as 

        	Area
        |--------->
        |
Area    |    x_{ij}  represents the number of common genes enriched in one structure over another
        |                *still not sure what this phrase means*
        V

'''

for column_label,column_data in df.iteritems():
	print column_label