# -*- coding: utf-8 -*-
# @Author  : zywang & ziye Wang
# @FileName: SolidBin.py


from sklearn.neighbors import kneighbors_graph
import pandas as pd
from scipy.sparse import csc_matrix,coo_matrix
import time
import numpy as np

from scipy.special import perm

from scipy.sparse import linalg
from sklearn.cluster import KMeans
from sklearn import metrics

import scipy.sparse as sp
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



'''def gen_knn_affinity_graph_large_scale(X, top=500):
    A = kneighbors_graph(X, n_neighbors=min(top, X.shape[0] - 1), mode='distance', p=1,
                         n_jobs=-1)
    nzIdx = np.nonzero(A)
    bandwidth = np.percentile(A[nzIdx[0], nzIdx[1]], 5)
    A[nzIdx[0], nzIdx[1]] = np.exp(-A[nzIdx[0], nzIdx[1]] / bandwidth)
    return A'''

def gen_knn_affinity_graph(X, top=np.Inf):
    A = kneighbors_graph(X, n_neighbors=min(top, X.shape[0] - 1), mode='distance', p=1,
                         n_jobs=-1)
    min_ = A.min()
    max_ = A.max()
    A=A.toarray()
    A = 1 - (A - min_)/(max_-min_)
    A =coo_matrix(A)
    return A


def gen_X(com_file, cov_file):
    covHeader = pd.read_csv(cov_file, sep='\t', nrows=1)
    covMat = pd.read_csv(cov_file, sep='\t', usecols=range(1, covHeader.shape[1])).as_matrix()
    namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]
    mapObj = dict(zip(namelist, range(len(namelist))))

    compositHeader = pd.read_csv(com_file, sep=',', nrows=1)
    shuffled_compositMat = pd.read_csv(com_file, sep=',', usecols=range(1, compositHeader.shape[1])).as_matrix()
    shuffled_namelist = pd.read_csv(com_file, sep=',', usecols=range(1)).values[:, 0]

    covIdxArr = []
    for contigIdx in range(len(shuffled_namelist)):
        assert (shuffled_namelist[contigIdx] in mapObj)
        covIdxArr.append(mapObj[shuffled_namelist[contigIdx]])
    compositMat = shuffled_compositMat.copy()
    compositMat[covIdxArr, :] = shuffled_compositMat

    covMat = covMat + 1e-2
    covMat = covMat / covMat.sum(axis=0)[None, :]
    if covMat.shape[1] >= 10:
        covMat = covMat / covMat.sum(axis=1)[:, None]
    compositMat = compositMat + 1
    compositMat = compositMat / compositMat.sum(axis=1)[:, None]
    X_t = np.hstack((covMat, compositMat)) * 1e1
    return X_t, namelist, mapObj


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

    if Dsqrt.ndim == 1:
        logger.error("matrix must have two dims")
        sys.exit(0)

    DW = sp.linalg.spsolve(Dsqrt, W)
    LT = sp.linalg.spsolve(Dsqrt.T, DW.T)
    L = (LT.T + LT) / 2
    d, v = linalg.eigs(L, k, which='LR')
    uu, dummy = np.real(v).T, np.real(d)
    return uu.astype(np.float)


def normalize_2nd(vv):
    n = vv.shape[1]
    ss = np.sqrt(np.sum(np.square(vv), axis=0))
    c = np.divide(vv, ss, out=np.zeros_like(vv), where=ss != 0)
    return c


def silhouette(X, W, label):
    X_colsum = np.sum(X ** 2, axis=1)
    X_colsum = X_colsum.reshape(len(X_colsum), 1)
    W_colsum = np.sum(W ** 2, axis=1)
    W_colsum = W_colsum.reshape(len(W_colsum), 1)

    Dsquare = np.tile(X_colsum, (1, W.shape[0])) + np.tile(W_colsum.T, (X.shape[0], 1)) - 2 * X.dot(W.T)
    # avoid error caused by accuracy
    Dsquare[Dsquare < 0] = 0
    D = np.sqrt(Dsquare)
    aArr = D[np.arange(D.shape[0]), label]
    D[np.arange(D.shape[0]), label] = np.inf
    bArr = np.min(D, axis=1)
    tmp = (bArr - aArr) / np.maximum(aArr, bArr)
    return np.mean(tmp)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def gen_bestk(contig_file, bestK=0):
    if bestK == 0:
        fragScanURL = os.path.join(os.getcwd(), 'auxiliary', 'FragGeneScan1.19', 'run_FragGeneScan.pl')
        os.system("chmod 777 " + fragScanURL)
        hmmExeURL = os.path.join(os.getcwd(), 'auxiliary', 'hmmer-3.1b1', 'bin', 'hmmsearch')
        os.system("chmod 777 " + hmmExeURL)
        markerExeURL = os.path.join(os.getcwd(), 'auxiliary', 'test_getmarker.pl')
        os.system("chmod 777 " + markerExeURL)
        markerURL = os.path.join(os.getcwd(), 'auxiliary', 'marker.hmm')
        seedURL = contig_file + ".seed"
        fragResultURL = contig_file + ".frag.faa"
        hmmResultURL = contig_file + ".hmmout"

        candK = 5
        maxK = 200
        stepK = 5

        if not (os.path.exists(fragResultURL)):
            fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + contig_file + ".frag -complete=0 -train=complete -thread=10 1>" + args.contig_file + ".frag.out 2>" + args.contig_file + ".frag.err"
            print("exec cmd: " + fragCmd)
            os.system(fragCmd)

        if os.path.exists(fragResultURL):
            if not (os.path.exists(hmmResultURL)):
                hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu 10 " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
                print("exec cmd: " + hmmCmd)
                os.system(hmmCmd)

            if os.path.exists(hmmResultURL):
                if not (os.path.exists(seedURL)):
                    markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " 1000 " + seedURL
                    print("exec cmd: " + markerCmd)
                    os.system(markerCmd)

                if os.path.exists(seedURL):
                    candK = file_len(seedURL)
                    maxK = 2 * candK
                    stepK = 2  
                else:
                    print("markerCmd failed! Not exist: " + markerCmd)
            else:
                print("Hmmsearch failed! Not exist: " + hmmResultURL)
        else:
            print("FragGeneScan failed! Not exist: " + fragResultURL)

        bestK = candK
        bestSilVal = 0
        t = time.time()
        for k in range(candK, maxK, stepK):
            kmeans = KMeans(n_clusters=k, init='k-means++')
            kmeans.fit(X_t)
            silVal = silhouette(X_t, kmeans.cluster_centers_, kmeans.labels_)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal) + "\telapsed time:" + str(time.time() - t))
            t = time.time()

            if silVal > bestSilVal:
                bestSilVal = silVal
                bestK = k
            else:
                break

        candK = bestK + 4
        bestSilVal_2nd = 0
        for k in range(candK, maxK, stepK):
            kmeans = KMeans(n_clusters=k, init='k-means++')
            kmeans.fit(X_t)
            silVal_2nd = silhouette(X_t, kmeans.cluster_centers_, kmeans.labels_)
            logger.info("k:" + str(k) + "\tsilhouette:" + str(silVal_2nd) + "\telapsed time:" + str(time.time() - t))
            t = time.time()
            if silVal_2nd > bestSilVal_2nd:
                bestSilVal_2nd = silVal_2nd
                bestK = k
            else:
                break
        if bestSilVal_2nd > bestSilVal:
            bestSilVal = bestSilVal_2nd
        else:
            bestK = candK - 4
        logger.info("bestk:" + str(bestK) + "\tsilVal:" + str(bestSilVal))

    else:
        logger.info("Use the pre-specified cluster number! k=" + str(bestK))

    return bestK


def getLinkRatio(label):
    cluster = np.unique(label)
    clusterNum = cluster.shape[0]
    linkNum = 0
    for i in range(clusterNum):
        ind = np.nonzero(label == i)
        linkNum = linkNum + perm(len(ind[0]), 2)
    ratio = linkNum / perm(label.shape[0], 2)
    return ratio


def buildSFSMat(A, mlRatio):
    simCutoff = np.percentile(A.toarray(), 100 - 100 * mlRatio,interpolation='midpoint')
    (a1, a2) = np.nonzero(A > simCutoff)
    ML = csc_matrix((np.ones(a1.shape[0]), (a1[:], a2[:])), shape=(A.shape[0], A.shape[0]))
    ML = maximum(ML, ML.T)
    ML = ML-np.diag(ML.diagonal())
    (a1, a2) = np.nonzero(ML)
    constraint_num=len(a1)/2
    logger.info("build SFS constraint num:\t"+str(constraint_num))
    return csc_matrix(ML)

def arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--contig_file', type=str,help=("The contigs file."))
    parser.add_argument('--coverage_profiles',type=str, help=(
        "The coverage profiles, containing a table where each row correspond to a contig, and each column correspond to a sample. All values are separated with tabs."))
    parser.add_argument('--composition_profiles',type=str, help=(
        "The composition profiles, containing a table where each row correspond to a contig, and each column correspond to the kmer composition of particular kmer. All values are separated with comma."))
    parser.add_argument('--priori_ml_list',type=str, help=(
        "The edges encoding either the co-alignment or other priori information, one row of an edge is in the format: contig_name_A\tcontig_name_B\t1. The edge is undirected."))
    parser.add_argument('--priori_cl_list',type=str, help=(
        "The cannot-link information, one row of an edge is in the format: contig_name_A\tcontig_name_B\t1. The edge is undirected."))
    parser.add_argument('--output',type=str, help="The output file, storing the binning result.")
    parser.add_argument('--log', type=str,help="Specify where to store log file")
    parser.add_argument('--clusters', default=0, type=int,
                        help="Specify the number of clusters. If not specified, the cluster number is estimated by single-copy genes.")
    parser.add_argument('-a', type=int, help="weight of must-link constraints")
    parser.add_argument('-b', type=int, help="weight of cannot-link constraints")
    parser.add_argument('--use_sfs',action="store_true",help="use SFS as constraints")


    args = parser.parse_args()

    if not (args.contig_file and args.coverage_profiles and args.composition_profiles and args.output):
        parser.error(
            "Data is missing, add file(s) using --contig_file <contig_file> and/or --coverage_profiles <abund_profiles> and/or --composition_profiles <comp_profiles> and/or --output <out_file>")
        sys.exit(0)
    return args

# ssf ssf-cl, use calinski_harabaz_score  (X)
# coalign cl, use calinski_harabaz_score (uut)
def ssncut_full_search(A, k, ml_mat, cl_mat, a_list, b_list,X=None):
    best_a = 0
    best_b = 0

    num_a = np.size(a_list)
    num_b = np.size(b_list)

    best_a_score=0
    best_b_score=0

    result = None

    if num_a > 1 or num_b==0:
        for alpha_idx in range(num_a):
            alpha = a_list[alpha_idx]
            labelPred, uu = ssncut(A, k, ml_mat, [], alpha, 0)
            if X is None:
                a_score = metrics.calinski_harabaz_score(uu.T, labelPred)
            else:
                a_score = metrics.calinski_harabaz_score(X, labelPred)
            logger.info("alpha:\t" + str(alpha) + "\tcalinski_harabaz_score_uut:\t" + str(a_score))

            if a_score > best_a_score:
                best_a_score=a_score
                result=labelPred
                best_a=alpha
        logger.info("estimated alpha is:\t" + str(best_a))

    if num_b > 0:

        alpha = best_a

        for beta_idx in range(num_b):
            beta = b_list[beta_idx]
            labelPred, uu = ssncut(A, k, ml_mat, cl_mat, alpha, beta)
            if X is None:
                b_score = metrics.calinski_harabaz_score(uu.T, labelPred)
            else:
                b_score = metrics.calinski_harabaz_score(X, labelPred)
            logger.info("alpha:\t" + str(alpha) + "\tbeta:\t" + str(beta) + "\tcalinski_harabaz_score_uut:\t" + str(
                b_score))

            if b_score>best_b_score:
                best_b_score=b_score
                result=labelPred
                best_b=beta

        logger.info("estimated beta is:\t" + str(best_b))

    return best_a, best_b,result


def ssncut(Simi_Matrix, K_Clust, ML, CL, alpha, beta):
    xx = mapping_all(Simi_Matrix, K_Clust, ML, CL, alpha, beta)
    kmeans = KMeans(n_clusters=K_Clust, init='k-means++')
    kmeans.fit(xx.T)
    return kmeans.labels_, xx

def save_result(result,filepath,namelist):
    filedir,filename=os.path.split(filepath)
    if not filename:
        filename="result.csv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    for contigIdx in range(len(namelist)):
        f.write(namelist[contigIdx] + "," + str(result[contigIdx].item(0)) + "\n")
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
    logger.info("priori_ml_list:\t" + (args.priori_ml_list if args.priori_ml_list else "Not Available"))
    logger.info("priori_cl_list:\t" + (args.priori_cl_list if args.priori_cl_list  else "Not Available"))
    logger.info("clusters:\t" + (str(args.clusters) if args.clusters > 0 else "Auto"))

    com_file = args.composition_profiles
    cov_file = args.coverage_profiles

    X_t, namelist, mapObj = gen_X(com_file, cov_file)
    contigNum = X_t.shape[0]
    start_time = time.time()
    affinity_mat = gen_knn_affinity_graph(X_t)
    affinity_mat = maximum(affinity_mat, affinity_mat.T)
    affinity_mat = csc_matrix(affinity_mat)
    end_time = time.time()
    logger.info("finish affinity matrix\ttime:"+str(end_time - start_time))
    contig_file = args.contig_file
    bestK = gen_bestk(contig_file, args.clusters)

    logger.info("estimate clusters:\t"+str(bestK))

    # search range for alpha and beta
    alpha_list = [1.0, 5.0, 10.0, 20.0, 30.0, 40.0] if not args.a else [args.a]
    sfs_alpha_list=[0]+alpha_list
    beta_list  = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0] if not args.b else [args.b]


    sfs_mat=None

    if args.use_sfs:
        logger.info("start building sfs matrix")
        oriLabels, uu = ssncut(affinity_mat, bestK, [], [], 0, 0)
        oriRatio = getLinkRatio(oriLabels)
        ratioUpper = (len(np.nonzero(affinity_mat)[0])) / perm(contigNum, 2)
        finalRatio = min(oriRatio * 0.4, ratioUpper)
        finalRatio = (finalRatio * perm(contigNum, 2) + contigNum) / (contigNum ** 2)
        logger.info("constraint link ratio:\t" + str(finalRatio))
        sfs_mat = buildSFSMat(affinity_mat, finalRatio)


    if not args.priori_ml_list and not args.priori_cl_list:
        if args.use_sfs:
            logger.info("SolidBin mode:\tSFS")
            a,b,res=ssncut_full_search(affinity_mat,bestK, sfs_mat,[] , sfs_alpha_list, [],X=X_t)
            logger.info("alpha:\t"+str(a))
            save_result(res,args.output,namelist)
        else:
            logger.info("SolidBin mode:\tnaive")
            res, uu = ssncut(affinity_mat, bestK, [], [], 0, 0)
            save_result(res, args.output, namelist)

    # currently, do not use both two constraints
    if args.priori_ml_list:
        if args.priori_cl_list:
            logger.info("two kinds of constraints provided, only use must-link constraints")
        mlList = pd.DataFrame(np.genfromtxt(args.priori_ml_list,
                                            converters={0: lambda s: mapObj[s.decode("utf-8")], 1: lambda s: mapObj[s.decode("utf-8")],
                                                        2: lambda s: float(s)}, delimiter=',')).as_matrix()
        coalignMat = csc_matrix((mlList[:, 2], (mlList[:, 0], mlList[:, 1])),
                                shape=(contigNum, contigNum))
        coalignMat = maximum(coalignMat, coalignMat.T)
        logger.info("SolidBin mode:\tcoalign")
        a,b,res=ssncut_full_search(affinity_mat, bestK, coalignMat, [], alpha_list, [])
        logger.info("alpha:\t" + str(a))
        save_result(res,args.output,namelist)

    if not args.priori_ml_list and args.priori_cl_list:
        clList = pd.DataFrame(np.genfromtxt(args.priori_cl_list, converters={0: lambda s: mapObj[s.decode("utf-8")],
                                                                                        1: lambda s: mapObj[s.decode("utf-8")],
                                                                                        2: lambda s: float(s)},
                                            delimiter=',')).as_matrix()
        cannotLinkMat = csc_matrix((clList[:, 2], (clList[:, 0], clList[:, 1])), shape=(contigNum, contigNum))
        cannotLinkMat = maximum(cannotLinkMat, cannotLinkMat.T)

        if args.use_sfs:
            logger.info("SolidBin mode\tSFS-CL")
            sfs_cl=sfs_mat-sfs_mat.multiply(cannotLinkMat)
            a,b,res=ssncut_full_search(affinity_mat, bestK, sfs_cl, [], sfs_alpha_list, [],X=X_t)
            logger.info("alpha:\t"+str(a))
            save_result(res,args.output,namelist)

        else:
            logger.info("SolidBin mode\tCL")
            a,b,res=ssncut_full_search(affinity_mat, bestK, [], cannotLinkMat, [], beta_list)
            logger.info("beta:\t"+str(b))
            save_result(res,args.output,namelist)


    if args.log:
        logger.removeHandler(handler)
    logger.removeHandler(console_hdr)
