#eg: ./run.sh <data_path>/input/strmgCAMI2_short_read_pooled_megahit_assembly.fasta 1000 4
# ./run.sh fasta_file  contig_length_threshold kmer_number


input_dir=$(dirname $1)

# prpocess.sh $1s
python scripts/filter_tooshort_for_contig_file.py $1 $2
python scripts/filter_tooshort_for_cov_file.py ${input_dir} $2
python scripts/gen_kmer.py $1 $2 $3
