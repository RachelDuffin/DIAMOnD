# DIAMOnD

Modification of DIAMOnD.py which runs on the output edge table from Cytoscape. Parses a cytoscape edge csv table and 
uses this as the input network_file to the DIAMOnD.py script. Alternatively the script can be run without cytoscape 
input, using a network PPI file as specified in the Example directory.

DIAMOnD.py runs the DIAMOnD algorithm as described in
 
 A DIseAse MOdule Detection (DIAMOnD) Algorithm derived from a
 systematic analysis of connectivity patterns of disease proteins in
 the Human Interactome. PlOS Comp Bio (in press), 2015.

by Susan Dina Ghiassian, Joerg Menche & Albert-Laszlo Barabasi

#--------------------
## Inputs
| Input                     | Details                                                                                                                                                                                                                                                                     |
|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| network_file              | Network file detailing interactions (edges) between nodes. This can either be an exported cytoscape level 0 edge table, or a delimiter-separated table. In the case of a delimiter-separated table, the first two columns are interpreted as interaction gene1 <==> gene2). |
| seed_file                 | Table containing the seed genes associated with a phenotype of interest, only column 1 is used.                                                                                                                                                                             |
| n                         | Number of diamond nodes to generate                                                                                                                                                                                                                                         |
| outfile_name              | Results will be saved under this file name                                                                                                                                                                                                                                  |
| alpha (optional)          | Supplies a weight to the seeds. Default value is 1.                                                                                                                                                                                                                         |
| cytoscape_file (optional) | True or False value. Denotes whether or not the analysis is being run using a cytoscape network file.                                                                                                                                                                       |

Note that gene IDs should be consistent between the network and seed files

# -------------------
## Running the code

Instruction to use the source code:
1. Download the code.
2. Make sure you are in the main directory where the code is.
3. Run the following.</br>
 <em><pre>python3 DIAMOnD.py  network_file seed_file  n  outfile_name  alpha(optional)  cytoscape_file(optional)</pre></em>

If you do not want to supply a weight to the seeds, or do not want to run the analysis using a cytoscape input file, 
leave alpha and cytoscape_file blank.

###Non-cytoscape mode
The Example directory contains two input files - seed_genes.txt and PPI.txt (protein-protein interaction network).

The following command will generate the first 100 DIAMOnD nodes and save them in a file:

<em><pre>python3  DIAMOnD.py  Example/PPI.txt  Example/seed_genes.txt  100  diamond_results.txt</pre></em>

###Cytoscape mode
The following command will generate the first 100 DIAMOnD nodes from a cytoscape input file and seed gene list and save
them in a file: 
<em><pre>python3  DIAMOnD.py  network_default_edge.csv  seed_genes.txt  100  diamond_results.txt True</pre></em>