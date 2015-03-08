#!/bin/sh

while read line; do
 	python ./src/make-figure-5-expansions.py -s "${line}" -o "./images/${line}-figure-5.tiff";
done < "$1"