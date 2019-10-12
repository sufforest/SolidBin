# -*- coding: utf-8 -*-
# @Author  :ziye Wang

import sys
import pandas as pd
from argparse import ArgumentParser
import numpy as np
from itertools import combinations

def arguments():
    parser = ArgumentParser()
    parser.add_argument('--TAXAassign_file', help=("The output file of TAXAassign"))

    args = parser.parse_args()
    if not (args.TAXAassign_file):
        parser.error(
            "Data is missing, add files(s) using --TAXAassign_file < The output file of TAXAassign > ")
        sys.exit(0)
    return args

def gen_coailgn_output(TAXAassign_file):
    namelist = pd.read_csv(TAXAassign_file, sep=',', header=None, usecols=range(1)).values[:, 0]
    taxaHeader = pd.read_csv(TAXAassign_file, sep=',', nrows=1)
    mapObj = dict(zip(namelist, range(len(namelist))))  #
    taxaassignMat = pd.read_csv(TAXAassign_file, sep=',', header=None, usecols=range(1, taxaHeader.shape[1])).values
    output_coalign_file = TAXAassign_file + '.coalign.csv'
    # label
    taxaassignMat_coalign = taxaassignMat.sum(1)
    taxa_label = np.unique(taxaassignMat_coalign)
    genomeNum = taxa_label.shape[0]
    contigs_cluster_map = [[] for i in range(genomeNum)]
    for i in range(genomeNum):
        for j in range(len(taxaassignMat_coalign)):  #
            if taxaassignMat_coalign[j] == taxa_label[i]:
                contigs_cluster_map[i].append(j)

    out_text = open(output_coalign_file, 'w')
    for i in range(genomeNum):
        temp = list(combinations(contigs_cluster_map[i], 2))
        for j in range(len(temp)):
            out_text.write(namelist[temp[j][0]] + ',' + namelist[temp[j][1]] + ',' + str(1))
            out_text.write('\n')

    out_text.close()

def gen_cannotlink_output(TAXAassign_file):
    namelist = pd.read_csv(TAXAassign_file, sep=',', header=None, usecols=range(1)).values[:, 0]
    taxaHeader = pd.read_csv(TAXAassign_file, sep=',', nrows=1)
    mapObj = dict(zip(namelist, range(len(namelist))))  #

    taxaassignMat = pd.read_csv(TAXAassign_file, sep=',', header=None, usecols=range(1, taxaHeader.shape[1])).values

    taxaassignMat_genus = taxaassignMat[:, 4]

    taxa_label_genus = np.unique(taxaassignMat_genus)
    genusNum = taxa_label_genus.shape[0]
    contigs_cluster_map_genus = [[] for i in range(genusNum)]

    output_cannot_file = TAXAassign_file + '.cannot.csv'
    out_text = open(output_cannot_file, 'w')
    for i in range(len(taxaassignMat_genus)):
        for j in range(i + 1, len(taxaassignMat_genus)):
            if taxaassignMat_genus[i] != taxaassignMat_genus[j]:
                out_text.write(namelist[i] + ',' + namelist[j] + ',' + str(1))
                out_text.write('\n')

    out_text.close()

if __name__ == '__main__':
    args = arguments()
    print("Input args:")
    print("TAXAassign_file:\t" + args.TAXAassign_file)
    print("Generate coalign file:")
    gen_coailgn_output(args.TAXAassign_file)
    print("Generate coalign file:\t Done!")
    print("Generate cannot link file:")
    gen_cannotlink_output(args.TAXAassign_file)
    print("Generate cannot link file:\t Done!")
    print("Finished!")

