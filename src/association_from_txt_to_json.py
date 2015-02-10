import json

from awesome_print import ap 


data = open('go_association.txt','rb').read().splitlines()

data_as_dict = {line.split()[0]:line.split()[-1].split(';') for line in data}




json.dump(data_as_dict,open('go_association.json','wb'))
