#!/usr/bin/python3
'''
'''
from __future__ import print_function
import numpy as np
from itertools import product
from Bio import SeqIO
# optimized sliding window function from
# http://stackoverflow.com/a/7636587
from itertools import tee
from collections import Counter, OrderedDict
import pandas as pd
import argparse
import mimetypes
import os

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in range(i):
            next(el, None)
    return zip(*els)

def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash,counter

def generate_features_from_fastx(fastx_file,length_threshold,kmer_len,outfile,length_file=None):
    kmer_dict,nr_features = generate_feature_mapping(kmer_len)
    file_type = mimetypes.guess_type(fastx_file)[1]
    if file_type == 'gzip':
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError("Unknown type of file: '{}".format(fastx_file))

    length = {}
    if os.path.getsize(fastx_file) == 0:
        return length

    file_format = None
    line = f.readline()
    if line.startswith('@'):
        file_format = "fastq"
    elif line.startswith(">"):
        file_format = "fasta"
    f.seek(0)

    if not file_format:
        raise RuntimeError("Invalid sequence file: '{}".format(fastx_file))
        
    # Store composition vectors in a dictionary before creating dataframe
    composition_d = OrderedDict()
    contig_lengths = OrderedDict()

    for seq in SeqIO.parse(fastx_file,file_format):
        seq_len = len(seq)
        if seq_len <= length_threshold:
            continue
        contig_lengths[seq.id] = seq_len
        # Create a list containing all kmers, translated to integers
        kmers = [
            kmer_dict[kmer_tuple]
            for kmer_tuple
            in window(str(seq.seq).upper(), kmer_len)
            if kmer_tuple in kmer_dict
        ]
        kmers.append(nr_features-1)
        composition_v = np.bincount(np.array(kmers,dtype=np.int64))
        composition_v[-1]-=1
        composition_d[seq.id] = composition_v 
    df = pd.DataFrame.from_dict(composition_d, orient='index', dtype=float)
    df.to_csv(outfile)
    if length_file:
        df = pd.DataFrame.from_dict(contig_lengths, orient='index', dtype=float)
        df.to_csv(length_file)

    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Extract k-mer and length")
    parser.add_argument("-f", "--fastx_file", help="FASTA or FASTQ file", required=True)
    parser.add_argument("-k", type=int,default=4)
    parser.add_argument("-t",type=int,default=1000,help="length threshold")
    parser.add_argument("-o",type=str,default="kmer.csv",help="output file")
    parser.add_argument("-l",type=str,default=None,help="length file")
    args = parser.parse_args()
    generate_features_from_fastx(args.fastx_file,args.t,args.k,args.o,args.l)
