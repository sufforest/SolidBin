# SolidBin
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
source activate cocacola_env
```

Install the SolidBin dependencies into this environment:

```sh
conda install numpy pandas scikit-learn scipy
```

## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate coverage profile and composition profile as input to our program.

There are several methods that can generate these two types of information.

We recommend using the way CONCOCT adopts, you can find it [here](https://concoct.readthedocs.io/en/latest/complete_example.html).

## <a name="usage"></a>Usage


> - Uage:          [--contig_file CONTIG_FILE]
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
                        





