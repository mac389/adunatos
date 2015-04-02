find ./ALE-output/ -name '*.structures' | while read fname;
 do 
 	python seaborn-figure-5.py -i $fname -c 100 
 done