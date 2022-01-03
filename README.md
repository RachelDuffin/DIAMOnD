# GenePanelAnalyzer

This script modifies DIAMOnD.py and makes use of the Phen2Gene tool. It carries out the following steps
1. Parses a cytoscape edge csv table for use as the input network file (optional)
2. Runs the DIAMOnD script on the network and seed file which considers the significance of the number of connections
to seed genes to prioritize the genes in the network for their disease relevance. The script has been modified to only 
output genes with a p value of <0.001
3. Parses the output of DIAMOnD and converts NCBI IDs to genenames using the Biomart python package
4. Inputs the seed and predicted genes, and the disease-associated HPO terms, to Phen2Gene to produce an output scored
candidate gene list.

#----------------------------------------------------------------------------------------------------------------------
## Inputs
| Input                     | Details                                                                                                                                                                                                                                                                     |
|---------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| hpo_terms                 | File containing HPO terms                                                                                                                                                                                                                                                   |
| network_file              | Network file detailing interactions (edges) between nodes. This can either be an exported cytoscape level 0 edge table, or a delimiter-separated table. In the case of a delimiter-separated table, the first two columns are interpreted as interaction gene1 <==> gene2). |
| seed_file                 | Table containing the seed gene names associated with a phenotype of interest, only column 1 is used.                                                                                                                                                                        |
| n                         | Number of diamond nodes to generate                                                                                                                                                                                                                                         |
| outfile_name              | Results of the DIAMOnD analysis will be saved under this file name                                                                                                                                                                                                          |
| alpha (optional)          | Supplies a weight to the seeds. Default value is 1.                                                                                                                                                                                                                         |
| cytoscape_file (optional) | True or False value. Denotes whether or not the analysis is being run using a cytoscape network file.                                                                                                                                                                       |

Note that gene IDs in the network file should be NCBI IDs.
#----------------------------------------------------------------------------------------------------------------------
## Installation and running the code
The Example directory contains a set of input files that can be used to run the code.

Instruction to use the source code:
1. Clone the repository, then install the requirements in a virtual environment.
 <em><pre>pip3 install -r requirements.txt</pre></em>
2. Install Phen2Gene by running the following and specifying the locations for install
 <em><pre>bash setup.sh</pre></em>
3. Make sure you are in the main root directory of the repository then run the commands specified under either of the
headings below

If you do not want to supply a weight to the seeds, or do not want to run the analysis using a cytoscape input file, 
leave alpha and cytoscape_file blank.

### Non-cytoscape mode
The following command will generate the first 100 DIAMOnD nodes and save them in a file:
<em><pre>python3 GenePanelAnalyzer.py Example/CM_core_HPO.txt Example/CM_PPI_file.txt Example/CM_green_genes_names.txt 100 CM_diamond_output.txt</pre></em>

### Cytoscape mode
The following command will generate the first 100 DIAMOnD nodes from a cytoscape input file and save them in a file: 
<em><pre>python3 GenePanelAnalyzer.py Example/CM_core_HPO.txt Example/CM_human_default_edge.csv Example/CM_green_genes_names.txt 100 CM_diamond_output.txt True</pre></em>

#----------------------------------------------------------------------------------------------------------------------
## Citations

[Zhao, M., Havrilla, J. M., Fang, L., Chen, Y., Peng, J., Liu, C., Wu C., Sarmady M., Botas P., Isla J., Lyon G., 
Weng C., Wang, K. (2019). Phen2Gene: Rapid Phenotype-Driven Gene Prioritization for Rare Diseases.NAR Genomics and 
Bioinformatics, Volume 2, Issue 2, June 2020, lqaa032](https://doi.org/10.1093/nargab/lqaa032)

[Ghiassian SD, Menche J, Barabási AL. A DIseAse MOdule Detection (DIAMOnD) algorithm derived from a systematic analysis 
of connectivity patterns of disease proteins in the human interactome. PLoS computational biology. 2015 Apr 
8;11(4):e1004120.](https://doi.org/10.1371/journal.pcbi.1004120.g001)
