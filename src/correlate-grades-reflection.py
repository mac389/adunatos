import csv

from awesome_print import ap 

data = list(csv.DictReader(open('data.csv','rb')))

ap(data)