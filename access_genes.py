import requests

from awesome_print import ap 
base = "http://api.brain-map.org/api/v2"
ontologyId = 1
structureGraphId = 1

sectionDataSetId = 69855739 #This seems like a magic constant?

#Download structure graph
structure_url = '/'.join([base,'structure_graph_download','%s.json'%structureGraphId])

structure_ontology = requests.get(structure_url).json()['msg']

ap(list(datum['ontology_id'] for datum in structure_ontology))


expression_url = '/'.join([base,"data/StructureUnionize","query.json"])
parameters = {"criteria":"[structures.ontology_id$eq%d][section_data_set_id$eq%d]"%(ontologyId,sectionDataSetId),
		"include":"structure","start_row":0,"num_rows":2000}


expression_values = requests.get(expression_url,params=parameters)
ap(expression_values.url)
ap(expression_values.json())
'''
	apiQuery("data/StructureUnionize/query.json" + 
			 "?criteria=[structures.ontology_id$eq" + ontologyId + "]" + 
			 "[section_data_set_id$eq" + id + "]" +
			 "&include=structure",
			 processExpression);


#requests for num_rows, total_rows
'''