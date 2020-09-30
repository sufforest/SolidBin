#eg: bash scripts/run.sh <data_path>/input/final.contigs.fa 1000 4
#bash scripts/run.sh fasta_file contig_length_threshold kmer_number


input_dir=$(dirname $1)

# prpocess.sh $1s
python scripts/filter_tooshort_for_contig_file.py $1 $2
python scripts/filter_tooshort_for_cov_file.py ${input_dir} $2
python scripts/gen_kmer.py $1 $2 $3
