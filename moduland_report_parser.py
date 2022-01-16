#!/bin/python3

import sys

def main():
    suid_file = sys.argv[1]
    moduland_report = sys.argv[2]
    out = sys.argv[3]

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
    out_file = []

    with open(suid_file) as suid_file:
        for line in suid_file:
            elements = line.split(',')
            for ID in suid_list:
                if ID == elements[0]:
                    out_file.append(elements[1] + "\n")

    with open(out, "w") as f:
        f.writelines(out_file)


if __name__ == "__main__":
    main()