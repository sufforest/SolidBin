# SolidBin: Improving Metagenome Binning with Semi-supervised Normalized Cut
A genome binning method for contig binning, based on semi-supervised spectral clustering method.

## <a name="started"></a>Getting Started

### <a name="docker"></a>Conda

We recommend using conda to run SolidBin. Download [here](https://www.continuum.io/downloads)

After installing Anaconda (or miniconda), fisrt obtain SolidBin:

```sh
git clone https://github.com/sufforest/SolidBin
```
Then simply create a solidbin environment 

```sh
conda env create -f environment.yml
source activate solidbin
```



### <a name="docker"></a>Docker

We also provide our docker image. If you are more familiar with docker, you can just get our image by:

```sh
docker pull sufforest/solidbin
```

A simple way to use our image is just mount you data directory and run:

```sh
docker run -it -v /data/StrainMock/input:/input -v /data/StrainMock/output:/output solidbin python SolidBin.py --contig_file /input/StrainMock_Contigs_cutup_10K_nodup_filter_1K.fasta --composition_profiles /input/kmer_4.csv --coverage_profiles /input/cov_inputtableR.tsv --output /output/result.tsv --log /output/log.txt
```

Suppose that /data/StrainMock contains the data in your machine, this command mount two directories into the container so that our SolidBin can use them.

If you do not have composition or coverage profiles,  you can just enter the container and generate them by yourself.

```sh
docker run -it -v /data/StrainMock/input:/input -v /data/StrainMock/output:/output solidbin sh
```

## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate coverage profile and composition profile as input to our program.

There are several methods that can generate these two types of information and we provide one of them below.


### Composition Profile

Composition profile is the vector representation of contigs and we use kmer to generate this information.

```
$ python scripts/gen_kmer.py /path/to/data/contig.fasta 4 /path/to/input/composition.csv
```
Here we choose k=4. By default we only keep contigs longer than 1000, you can specify a different number using -t.

### Coverage Profile
For the coverage profile we use minimap, since it can address both short read and long read samples.

If you use SolidBin in docker, all dependencies and environment variables have been configured correctly, you can just mount you input directory into /input in docker, then slightly modify gen_cov.sh and run it.

You input directory should look like this:

```
.
+-- assembly.fasta
+-- sr
|   +-- short_read_sample_1
|   +-- short_read_sample_2
+-- pb
|   +-- pacbio_sample_1
|   +-- pacbio_sample_2
|   +-- pacbio_sample_3
```

For conda environment , you should check whether perl is installed.



## <a name="usage"></a>Usage


> - Usage:         [--contig_file CONTIG_FILE]
                   [--coverage_profiles COVERAGE_PROFILES]
                   [--composition_profiles COMPOSITION_PROFILES]
                   [--priori_ml_list ML_LIST] 
                   [--priori_cl_list CL_LIST] 
                   [--output OUTPUT]
                   [--log LOG_LOCATION]
                   [--clusters CLUSTERS]
                   [-a ALPHA]
                   [-b BETA]
                   [--use_sfs]

> - arguments

  	--contig_file CONTIG_FILE: 
              The contigs file.
	
  	--coverage_profiles COVERAGE_PROFILES: 
              The coverage_profiles, containing a table where each
              row correspond to a contig, and each column correspond
              to a sample. All values are separated with tabs.
  	--composition_profiles COMPOSITION_PROFILES: 
              The composition profiles, containing a table where
              each row correspond to a contig, and each column
              correspond to the kmer composition of particular kmer.
              All values are separated with comma.
	
  	--output OUTPUT:
              The output file, storing the binning result.

> - optional

  	--priori_ml_list ML_LIST:
              The list of contig pairs under must-link constraints, one row for one constraint in
              the format: contig_1,contig_2,weight.
                        
  	--priori_cl_list CL_LIST:
              The list of contig pairs under cannot-link constraints, one row for one constraint in
              the format: contig_1,contig_2,weight.
    --log LOG_LOCATION:
              The log file, storing the binning process and parameters.
    
    --clusters CLUSTERS: 
              Specify the number of clusters. If not specified, the
              cluster number is estimated by single-copy genes.
                        
    -a ALPHA:
              Specify the parameter of must-link constraints, if not specified, 
              alpha is estimated by an one-dimensional search.
    -b BETA:
              Specify the parameter of cannot-link constraints, if not specified, 
              beta is estimated by an one-dimensional search.
                  
    --use_sfs:
              Use sequence feature similarity to generate must-link constraints.
> - different SolidBin modes

  SolidBin mode | Usage  
  ------------- | -------------
 SolidBin-naive | --contig_file --coverage_profiles --composition_profiles --output --log 
 SolidBin-SFS   | --contig_file --coverage_profiles --composition_profiles --output --log --use_sfs
 SolidBin-coalign   | --contig_file --coverage_profiles --composition_profiles --output --log --priori_ml_list
 SolidBin-CL   | --contig_file --coverage_profiles --composition_profiles --output --log --priori_cl_list
 SolidBin-SFS-CL   | --contig_file --coverage_profiles --composition_profiles --output --log --priori_cl_list --use_sfs

## <a name="preprocessing"></a>MATLAB VERSION
We also provide MATLAB version code for the reproduction of the results in our paper "SolidBin: Improving Metagenome Binning with Semi-supervised Normalized Cut".

We have wrapped our MATLAB code in Python and you can use it in the same way as Python (Please copy "auxiliary" into "matlab_ver" before you use MATLAB-ver SolidBin).

NOTICE: We only test this wrapper on MATLAB R2017a and we can not guarantee the stabibility on other versions of MATLAB.

## <a name="preprocessing"></a>Contacts and bug reports
Please send bug reports or questions (such as the appropriate modes for your datasets) to
Ziye Wang: zwang17@fudan.edu.cn and Dr. Shanfeng Zhu: zhusf@fudan.edu.cn

## <a name="preprocessing"></a>References

[1] Lu, Yang Young, et al. "COCACOLA: binning metagenomic contigs using sequence COmposition, read CoverAge, CO-alignment and paired-end read LinkAge." Bioinformatics 33.6 (2017): 791-798.

[2] Alneberg, Johannes, et al. "Binning metagenomic contigs by coverage and composition." Nature methods 11.11 (2014): 1144.             

[3] Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome Research, 25: 1043–1055.

[4] Graham ED, Heidelberg JF, Tully BJ. (2017) "BinSanity: unsupervised clustering of environmental microbial assemblies using coverage and affinity propagation." PeerJ 5:e3035

## <a name="preprocessing"></a>Citation
Wang Z., et al. "SolidBin: Improving Metagenome Binning with Semi-supervised Normalized Cut." Bioinformatics. 2019 Apr 12. pii: btz253. doi: 10.1093/bioinformatics/btz253.

Note: If you cite SolidBin in your paper, please specify the mode of SolidBin you used in your paper to avoid confusion. Thank you. 



