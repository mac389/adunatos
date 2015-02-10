# adunatos


##### Quickstart
1. To create gene-structure matrix from Allen Brain Atlas

     python access-aba.py

This module saves its output to as JSON and TXT. The module writes to the TXT file after each query. It dumps the JSON file after it has finished querying the Allen Brain Atlas's servers. 

To run
	
     cd adunatos
     ./main.sh

If main not recognized as executable replace the last line with these two lines:

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
