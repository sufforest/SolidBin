# SolidBin: Improving Metagenome Binning with Semi-supervised Normalized Cut
A genome binning method for contig binning, based on semi-supervised spectral clustering method.

## <a name="started"></a>Getting Started

We recommend using Anaconda to run SolidBin. Download [here](https://www.continuum.io/downloads)

After installing Anaconda, fisrt obtain SolidBin:

```sh
git clone https://github.com/sufforest/SolidBin
```
Create a Anaconda environment and activate it:

```sh
conda create -n solidbin python=3.6
source activate solidbin
```

Install the SolidBin dependencies into this environment:

```sh
$ conda install numpy pandas scikit-learn scipy
```

## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate coverage profile and composition profile as input to our program.

There are several methods that can generate these two types of information.

We recommend using the way CONCOCT adopts, you can find it [here](https://github.com/BinPro/CONCOCT/). 

For the reproducibility, we provide a complete example to show how each profile is generated:

### Dependency

Fisrt of all, we download [CONCOCT v4.0](https://github.com/BinPro/CONCOCT/archive/0.4.0.zip) and set environmental variables for software dependencies.

```sh
$ CONCOCT=/path/to/your/CONCOCT
$ MRKDUP=/path/to/your/picard-tools-1.77/MarkDuplicates.jar

```

We use Bowtie2 to map the reads of each sample back to the assembly.
And then we use MarkDuplicates to remove PCR duplicates, use BEDTools to compute the coverage histogram for each bam file.

After we install all dependecies and set environmental variables correctly, we can generate coverage profile and composition profile for the given dataset(Make sure that directories of dependency are both in the PATH environmental variable).

Take dataset *StrainMock* as an example.


### Composition Profile

Composition profile is the vector representation of contigs and we use kmer(k=4) to generate this information.

```
$ python $CONCOCT/scripts/fasta_to_features.py /path/to/data/StrainMock_Contigs_cutup_10K_nodup_filter_1K.fasta 9417 4 /path/to/input/composition.csv
```
9417 is the number of contigs in this data and 4 is the k we choose.

### Coverage Profile
For the coverage profile, we first create the index on the assembly for bowtie2.

```
$ bowtie2-build StrainMock_Contigs_cutup_10K_nodup_filter_1K.fasta StrainMock_Contigs_cutup_10K_nodup_filter_1K.fasta
```
Then we map the reads of each sample (Here we choose *Sample1006*).

```
$CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 10 -p '-f' /path/to/samples/Sample1006_1.fasta /path/to/samples/Sample1006/Sample1006_2.fasta pair /path/to/data/StrainMock_Contigs_cutup_10K_nodup_filter_1K.fasta Sample1006 /path/to/samples/Sample1006/
```

Finally, we can generate coverage profile from the bam files for each sample.

```
python $CONCOCT/scripts/gen_input_table.py \
/path/to/data/StrainMock_Contigs_cutup_10K_nodup_filter_1K.fast \
/path/to/samples/Sample1006/Sample1006_pair-smds.bam /path/to/other/sample/bamfile... \
> /path/to/input/coverage.tsv
```



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

## <a name="preprocessing"></a>Citation
Wang Z., et al. "SolidBin: Improving Metagenome Binning with Semi-supervised Normalized Cut." Bioinformatics. 2019 Apr 12. pii: btz253. doi: 10.1093/bioinformatics/btz253.

Note: If you cite SolidBin in your paper, please specify the mode of SolidBin you used in your paper to avoid confusion. Thank you. 



