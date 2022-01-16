#!/bin/python3

import sys
from GenePanelAnalyzer import ncbi_id_to_genename, file_to_list

def main():
    suid_file = sys.argv[1]
    moduland_report = sys.argv[2]
    gene_name_file = sys.argv[3]
    out = sys.argv[4]

    # for line after
    suid_list = []

    with open(moduland_report) as f:
        for line in f:
            if line.startswith("the 10 core nodes"):
                nextline = next(f, '').strip('\n').strip(',')
                ids = nextline.split(',')
                for id in ids:
                    suid_list.append(id)

    suid_list = set(suid_list)
    final_id_list = []

    with open(suid_file) as suid_file:
        for line in suid_file:
            elements = line.split(',')
            for ID in suid_list:
                if ID == elements[0]:
                    final_id_list.append(elements[1])
    print(final_id_list)
    predicted_gene_names = ncbi_id_to_genename(final_id_list, "predicted gene names")
    print("PRed")
    print(predicted_gene_names)
    green_gene_namelist = file_to_list(gene_name_file)
    print("green")
    print(green_gene_namelist)
    combined_genelists = green_gene_namelist + predicted_gene_names
    print("combined")
    print(combined_genelists)

    with open(out, "w") as file:
        for gene in combined_genelists:
            file.write(gene + "\n")

if __name__ == "__main__":
    main()