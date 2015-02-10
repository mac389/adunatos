# adunatos


#### Quickstart
##### To _create_ gene-structure matrix from Allen Brain Atlas

     python access-aba.py

This module saves its output to as JSON and TXT. The module writes to the TXT file after each query. It dumps the JSON file after it has finished querying the Allen Brain Atlas's servers. 

##### To _isolate_ and _compare_ regions of interest
	
     cd adunatos
     ./main.sh

This shell script calls four python files to:
1. <b>divide</b> the gene-structure matrix into a control group that contains all brain areas and a regions of interest group (ROI), that contains all brain areas in which we are looking for enrichment. The file _filter_words_ specifies which structures to include in the regions of interest group. The output is a Pandas DataFrame.
1. <b>Annotate</b> the DataFrame with the Uniprot IDs. The DataFrame initially only has Entrez gene names. 
1. <b>Annotate</b> the DataFrame 
If your local machine does not recognize _main.sh_ as executable, replace the last line with these two lines:

     chmod +x ./main.sh
     ./main.sh

(Built and tested in Python 2.7)
#####Dependencies

- NumPy (built with 1.9.1)
- SciPy (built with 0.14.0)
- Matplotlib (1.4.0)
- Pandas (built with 0.14.1)
- mygene (built with 2.2.0)
- Awesome_print (0.1.3)
