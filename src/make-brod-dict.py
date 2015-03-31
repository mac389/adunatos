import json

from awesome_print import ap 
txt=open('broadman->structure.txt','rb').read().splitlines()

jdict = dict([line.split(':') for line in txt])
json.dump(jdict,open('broadman.json','wb'))