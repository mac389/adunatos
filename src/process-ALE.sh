for name in ./ALE-output/*.txt
do 
	python 'ALE->list-of-areas.py' -i $name > tmp
done