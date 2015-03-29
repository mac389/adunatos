import utils as tech
import pandas as pd 
import json

df = pd.read_json('../docs/gene-expression-by-area2.json').dropna(axis=1)
gene_names = df.index.values
functions = {gene_name:tech.get_gene_function(gene_name,levels='all') for gene_name in gene_names}
json.dump(functions,open('../docs/gene-function.json','wb'))