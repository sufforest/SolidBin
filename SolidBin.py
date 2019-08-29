# -*- coding: utf-8 -*-
# @Author  : zywang & ziye Wang
# @FileName: SolidBin.py


from sklearn.neighbors import kneighbors_graph
from sklearn.cluster import KMeans
from sklearn import metrics
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy import linalg
from scipy.special import perm
import scipy.sparse as sp
import mimetypes
from Bio import SeqIO
import logging
import os
import argparse
import sys

logger = logging.getLogger('SolidBin')

logger.setLevel(logging.INFO)

# logging
formatter = logging.Formatter('%(asctime)s - %(message)s')

console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)

logger.addHandler(console_hdr)


def gen_knn_affinity_graph(X, top=np.Inf):
    A = kneighbors_graph(X, n_neighbors=min(top, X.shape[0] - 1), mode='distance', p=1,
                         n_jobs=-1)
    min_ = A.min()
    max_ = A.max()
    A = A.toarray()
    A = 1 - (A - min_)/(max_-min_)
    A = csc_matrix(A)
    return A


def gen_X(com_file, cov_file):

    cov = pd.read_csv(cov_file, sep='\t', index_col=0)
    comp = pd.read_csv(com_file, sep=',', index_col=0)

    cov = cov.apply(lambda x: x+1e-5)
    cov = cov.div(cov.sum(axis=0), axis=1)
    cov = cov.div(cov.sum(axis=1), axis=0)

    comp = comp.apply(lambda x: x+1e-5)
    comp = comp.div(comp.sum(axis=1), axis=0)

    data = pd.merge(comp, cov, how='inner', on=None, left_on=None, right_on=None,
                    left_index=True, right_index=True, sort=True,
                    suffixes=('_comp', '_cov'), copy=True)

    namelist = data.index.tolist()
    mapObj = dict(zip(namelist, range(len(namelist))))

    return data.values, namelist, mapObj


def maximum(A, B):
    BisBigger = A - B
    BisBigger.data = np.where(BisBigger.data < 0, 1, 0)
    return A - A.multiply(BisBigger) + B.multiply(BisBigger)


def mapping_all(S, k, ML, CL, alpha, beta):
    D = sp.diags(np.asarray(S.sum(0))[0], format="csc")
    Dsqrt = D.sqrt()
    W = S
    if not isinstance(ML, list) and alpha:
        ML = csc_matrix(ML)
        alpha = alpha * sp.eye(S.shape[0])
        W = W - alpha * (sp.diags(np.asarray(ML.sum(0))[0]) - ML)

    if not isinstance(CL, list) and beta:
        CL = csc_matrix(CL)
        W = W - beta * CL

    DW = linalg.solve(Dsqrt.toarray(), W.toarray())
    LT = linalg.solve(Dsqrt.T.toarray(), DW.T)
    L = (LT.T + LT) / 2
    numL = L.shape[0]
    d, v = linalg.eigh(L, eigvals=(numL-k, numL-1))
    uu, dummy = np.real(v), np.real(d)
    return uu.astype(np.float)


def get_seed(contig_file, hard=0):
    if os.path.exists(contig_file+'.seed'):
        seed_list = []
        with open(contig_file+'.seed') as f:
            for line in f:
                seed_list.append(line.rstrip('\n'))
        return seed_list

    dir_path = os.path.dirname(os.path.realpath(__file__))
    fragScanURL = os.path.join(
        dir_path, 'auxiliary', 'FragGeneScan1.19', 'run_FragGeneScan.pl')
    os.system("chmod 777 " + fragScanURL)
    hmmExeURL = os.path.join(dir_path, 'auxiliary',
                             'hmmer-3.1b1', 'bin', 'hmmsearch')
    os.system("chmod 777 " + hmmExeURL)
    markerExeURL = os.path.join(dir_path, 'auxiliary', 'test_getmarker.pl')
    os.system("chmod 777 " + markerExeURL)
    markerURL = os.path.join(dir_path, 'auxiliary', 'marker.hmm')
    seedURL = contig_file+".seed"
    fragResultURL = contig_file+".frag.faa"
    hmmResultURL = contig_file+".hmmout"
    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL+" -genome="+contig_file+" -out="+contig_file + \
            ".frag -complete=0 -train=complete -thread=10 1>" + \
            contig_file+".frag.out 2>"+contig_file+".frag.err"
        print("exec cmd: "+fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL+" --domtblout "+hmmResultURL+" --cut_tc --cpu 10 " + \
                markerURL+" "+fragResultURL+" 1>"+hmmResultURL+".out 2>"+hmmResultURL+".err"
            print("exec cmd: "+hmmCmd)
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL+" "+hmmResultURL+" " + \
                    contig_file+" 1000 "+seedURL+" "+str(hard)
                print("exec cmd: "+markerCmd)
                os.system(markerCmd)

            if os.path.exists(seedURL):
                seed_list = []
                with open(seedURL) as f:
                    for line in f:
                        seed_list.append(line.rstrip('\n'))
                return seed_list
            else:
                return None
        else:
            print("Hmmsearch failed! Not exist: "+hmmResultURL)
    else:
        print("FragGeneScan failed! Not exist: "+fragResultURL)


def getLinkRatio(label):
    cluster = np.unique(label)
    clusterNum = cluster.shape[0]
    linkNum = 0
    for i in range(clusterNum):
        ind = np.nonzero(label == i)
        linkNum = linkNum + perm(len(ind[0]), 2)
    ratio = linkNum / perm(label.shape[0], 2)
    return ratio

def get_length(fastx_file):
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

    for seq_record in SeqIO.parse(f, file_format):
        length[seq_record.id] = len(seq_record.seq)

    f.close()

    return length

def buildSFSMat(A, mlRatio):
    simCutoff = np.percentile(
        A.toarray(), 100 - 100 * mlRatio, interpolation='midpoint')
    (a1, a2) = np.nonzero(A > simCutoff)
    ML = csc_matrix((np.ones(a1.shape[0]), (a1[:], a2[:])), shape=(
        A.shape[0], A.shape[0]))
    ML = maximum(ML, ML.T)
    ML = ML-np.diag(ML.diagonal())
    (a1, a2) = np.nonzero(ML)
    constraint_num = len(a1)/2
    logger.info("build SFS constraint num:\t"+str(constraint_num))
    return csc_matrix(ML)


def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--contig_file', type=str, help=("The contigs file."))
    parser.add_argument('--coverage_profiles', type=str, help=(
        "The coverage profiles, containing a table where each row correspond to a contig, and each column correspond to a sample. All values are separated with tabs."))
    parser.add_argument('--composition_profiles', type=str, help=(
        "The composition profiles, containing a table where each row correspond to a contig, and each column correspond to the kmer composition of particular kmer. All values are separated with comma."))
    parser.add_argument('--priori_ml_list', type=str, help=(
        "The edges encoding either the co-alignment or other priori information, one row of an edge is in the format: contig_name_A\tcontig_name_B\t1. The edge is undirected."))
    parser.add_argument('--priori_cl_list', type=str, help=(
        "The cannot-link information, one row of an edge is in the format: contig_name_A\tcontig_name_B\t1. The edge is undirected."))
    parser.add_argument('--output', type=str,
                        help="The output file, storing the binning result.(or a directory endswith '/') ")
    parser.add_argument('--log', type=str,
                        help="Specify where to store log file")
    parser.add_argument('--clusters', default=0, type=int,
                        help="Specify the number of clusters. If not specified, the cluster number is estimated by single-copy genes. By default, we use single-copy genes to initialize binning centers")
    parser.add_argument('-a', type=int, help="weight of must-link constraints")
    parser.add_argument(
        '-b', type=int, help="weight of cannot-link constraints")
    parser.add_argument('--use_sfs', action="store_true",
                        help="use SFS as constraints")

    args = parser.parse_args()

    if not (args.contig_file and args.coverage_profiles and args.composition_profiles and args.output):
        parser.error(
            "Data is missing, add file(s) using --contig_file <contig_file> and/or --coverage_profiles <abund_profiles> and/or --composition_profiles <comp_profiles> and/or --output <out_file>")
        sys.exit(0)
    return args


def ssncut_full_search(A, k, ml_mat, cl_mat, a_list, b_list, X=None):
    best_a = 0
    best_b = 0

    num_a = np.size(a_list)
    num_b = np.size(b_list)

    best_a_score = 0
    best_b_score = 0

    result = None

    if num_a > 1 or num_b == 0:
        for alpha_idx in range(num_a):
            alpha = a_list[alpha_idx]
            labelPred, uu = ssncut(A, k, ml_mat, [], alpha, 0)
            if X is None:
                a_score = metrics.calinski_harabasz_score(uu, labelPred)
            else:
                a_score = metrics.calinski_harabasz_score(X, labelPred)
            logger.info("alpha:\t" + str(alpha) +
                        "\tcalinski_harabasz_score_uut:\t" + str(a_score))

            if a_score > best_a_score:
                best_a_score = a_score
                result = labelPred
                best_a = alpha
        logger.info("estimated alpha is:\t" + str(best_a))

    if num_b > 0:

        alpha = best_a

        for beta_idx in range(num_b):
            beta = b_list[beta_idx]
            labelPred, uu = ssncut(A, k, ml_mat, cl_mat, alpha, beta)
            if X is None:
                b_score = metrics.calinski_harabasz_score(uu, labelPred)
            else:
                b_score = metrics.calinski_harabasz_score(X, labelPred)
            logger.info("alpha:\t" + str(alpha) + "\tbeta:\t" + str(beta) + "\tcalinski_harabasz_score_uut:\t" + str(
                b_score))

            if b_score > best_b_score:
                best_b_score = b_score
                result = labelPred
                best_b = beta

        logger.info("estimated beta is:\t" + str(best_b))

    return best_a, best_b, result


def ssncut(Simi_Matrix, K_Clust, ML, CL, alpha, beta):
    xx = mapping_all(Simi_Matrix, K_Clust, ML, CL, alpha, beta)
    if seed_list:
        assert len(seed_list) == K_Clust
        seed_idx=[mapObj[seed_name] for seed_name in seed_list]
        seed_center=xx[seed_idx][:]
        kmeans=KMeans(n_clusters=K_Clust, init=seed_center,n_init=1)
    else:
        kmeans = KMeans(n_clusters=K_Clust, init='k-means++')
    weight=[ length_dict[name] for name in namelist]
    weight=np.array(weight)
    weight=weight.astype(np.float32)
    kmeans.fit(xx,sample_weight=weight)
    return kmeans.labels_, xx


def save_result(result, filepath, namelist):
    filedir, filename = os.path.split(filepath)
    if not filename:
        filename = "result.csv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    f.write("@Version:0.9.0\n")
    f.write("@SampleID:{}\n".format(os.path.basename(contig_file)))
    f.write("@@SEQUENCEID\tBINID\n")
    for contigIdx in range(len(namelist)):
        f.write("{}\t{}\n".format(namelist[contigIdx],result[contigIdx].item(0)))
    f.close()


if __name__ == '__main__':
    args = arguments()

    if args.log:
        handler = logging.FileHandler(args.log)
        handler.setLevel(logging.INFO)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    logger.info("Input arguments:")
    logger.info("contig_file:\t" + args.contig_file)
    logger.info("coverage_profiles:\t" + args.coverage_profiles)
    logger.info("composition_profiles:\t" + args.composition_profiles)
    logger.info("output path:\t" + args.output)
    logger.info("priori_ml_list:\t" +
                (args.priori_ml_list if args.priori_ml_list else "Not Available"))
    logger.info("priori_cl_list:\t" +
                (args.priori_cl_list if args.priori_cl_list else "Not Available"))
    logger.info("clusters:\t" + (str(args.clusters)
                                 if args.clusters > 0 else "Auto"))

    com_file = args.composition_profiles
    cov_file = args.coverage_profiles

    X_t, namelist, mapObj = gen_X(com_file, cov_file)
    contigNum = X_t.shape[0]
    logger.info("Start generating affinity matrix")
    affinity_mat = gen_knn_affinity_graph(X_t)
    affinity_mat = maximum(affinity_mat, affinity_mat.T)
    affinity_mat = csc_matrix(affinity_mat)
    logger.info("Finish generating affinity matrix")
    contig_file = args.contig_file
    length_dict =get_length(contig_file)    
    seed_list=None
    if args.clusters==0:
        seed_list=get_seed(contig_file)
        bestK=len(seed_list)
    else:
        bestK = args.clusters

    logger.info("estimate clusters:\t"+str(bestK))

    # search range for alpha and beta
    alpha_list = [1.0, 5.0, 10.0, 20.0,40.0] if not args.a else [args.a]
    sfs_alpha_list = [0]+alpha_list
    beta_list = [0.02, 0.05, 0.1, 0.2] if not args.b else [args.b]

    sfs_mat = None

    if args.use_sfs:
        logger.info("start building sfs matix")
        oriLabels, uu = ssncut(affinity_mat, bestK, [], [], 0, 0)
        oriRatio = getLinkRatio(oriLabels)
        ratioUpper = (len(np.nonzero(affinity_mat)[0])) / perm(contigNum, 2)
        finalRatio = min(oriRatio * 0.4, ratioUpper)
        finalRatio = (finalRatio * perm(contigNum, 2) +
                      contigNum) / (contigNum ** 2)
        logger.info("constraint link ratio:\t" + str(finalRatio))
        sfs_mat = buildSFSMat(affinity_mat, finalRatio)

    if not args.priori_ml_list and not args.priori_cl_list:
        if args.use_sfs:
            logger.info("SolidBin mode:\tSFS")
            a, b, res = ssncut_full_search(
                affinity_mat, bestK, sfs_mat, [], sfs_alpha_list, [], X=X_t)
            logger.info("alpha:\t"+str(a))
            save_result(res, args.output, namelist)
        else:
            logger.info("SolidBin mode:\tnaive")
            res, uu = ssncut(affinity_mat, bestK, [], [], 0, 0)
            save_result(res, args.output, namelist)

    # currently, we do not use both two constraints
    if args.priori_ml_list:
        if args.priori_cl_list:
            logger.info(
                "two kinds of constraints provided, only use must-link constraints")
        mlList = pd.DataFrame(np.genfromtxt(args.priori_ml_list,
                                            converters={0: lambda s: mapObj[s.decode("utf-8")], 1: lambda s: mapObj[s.decode("utf-8")],
                                                        2: lambda s: float(s)}, delimiter=',')).as_matrix()
        coalignMat = csc_matrix((mlList[:, 2], (mlList[:, 0], mlList[:, 1])),
                                shape=(contigNum, contigNum))
        coalignMat = maximum(coalignMat, coalignMat.T)
        logger.info("SolidBin mode:\tcoalign")
        a, b, res = ssncut_full_search(
            affinity_mat, bestK, coalignMat, [], alpha_list, [])
        logger.info("alpha:\t" + str(a))
        save_result(res, args.output, namelist)

    if not args.priori_ml_list and args.priori_cl_list:
        clList = pd.DataFrame(np.genfromtxt(args.priori_cl_list, converters={0: lambda s: mapObj[s.decode("utf-8")],
                                                                             1: lambda s: mapObj[s.decode("utf-8")],
                                                                             2: lambda s: float(s)},
                                            delimiter=',')).as_matrix()
        cannotLinkMat = csc_matrix(
            (clList[:, 2], (clList[:, 0], clList[:, 1])), shape=(contigNum, contigNum))
        cannotLinkMat = maximum(cannotLinkMat, cannotLinkMat.T)

        if args.use_sfs:
            logger.info("SolidBin mode\tSFS-CL")
            sfs_cl = sfs_mat-sfs_mat.multiply(cannotLinkMat)
            a, b, res = ssncut_full_search(
                affinity_mat, bestK, sfs_cl, [], sfs_alpha_list, [], X=X_t)
            logger.info("alpha:\t"+str(a))
            save_result(res, args.output, namelist)

        else:
            logger.info("SolidBin mode\tCL")
            a, b, res = ssncut_full_search(
                affinity_mat, bestK, [], cannotLinkMat, [], beta_list)
            logger.info("beta:\t"+str(b))
            save_result(res, args.output, namelist)

    if args.log:
        logger.removeHandler(handler)
    logger.removeHandler(console_hdr)
