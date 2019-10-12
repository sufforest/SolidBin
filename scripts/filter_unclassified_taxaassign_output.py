import sys
import pandas as pd
from argparse import ArgumentParser
import numpy as np
def arguments():
    parser = ArgumentParser()
    parser.add_argument('--TAXAassign_file', help=("The output file of TAXAassign"))

    args = parser.parse_args()
    if not (args.TAXAassign_file):
        parser.error(
            "Data is missing, add files(s) using --TAXAassign_file < The output file of TAXAassign > ")
        sys.exit(0)
    return args
if __name__ == '__main__':
    args = arguments()
    print("Input args:")
    print("TAXAassign_file:\t" + args.TAXAassign_file)
    #TAXAassign_file = '/mnt/data4/wzy/preprocess/strmg_megahit_assembly/test_input/strmgCAMI2_short_read_pooled_gold_standard_assembly_ASSIGNMENTS.csv'
    TAXAassign_file=args.TAXAassign_file
    namelist = pd.read_csv(TAXAassign_file, sep=',', header=None, usecols=range(1)).values[:, 0]
    taxaHeader = pd.read_csv(TAXAassign_file, sep=',', nrows=1)
    mapObj = dict(zip(namelist, range(len(namelist))))  #
    taxaassignMat = pd.read_csv(TAXAassign_file, sep=',', header=None, usecols=range(0, taxaHeader.shape[1]))
    ex_list = list(taxaassignMat[6])
    ex_list_new = []
    for i in ex_list:
        if i != '__Unclassified__':
            ex_list_new.append(i)
    taxaassignMat = taxaassignMat[taxaassignMat[6].isin(ex_list_new)]
    output_file = TAXAassign_file + '.filter_unclassified.csv'
    taxaassignMat.to_csv(output_file, sep=',', index=None, header=None)


