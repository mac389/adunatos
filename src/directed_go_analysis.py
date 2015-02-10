import mygene, json

#GO for (DA, GABA, Glu, Inflammation, HAC/HDAC)
from awesome_print import ap 
from bioservices.kegg import KEGG 

#These number from Amigo
pathways_of_interest = {'DA':'GO:0007212','GABA':'GO:0007214',
  'Glu':'GO:0007215','Inflammation':"GO:0006954",
  'Epigenetics':"GO:0040029"}

#get proteins for each paths

mg = mygene.MyGeneInfo()

uniprot_go = json.load(open('uniprot->goa.json','rb'))

for item in uniprot_go:
	ap(item)

#for moniker,go in pathways_of_interest.iteritems():
