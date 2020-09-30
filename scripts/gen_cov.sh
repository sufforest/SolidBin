
# An example to generate coverage for hybrid dataset

# input_dir: the directory containts all data we need
# assembly: fasta file consisting of contigs
# mapdir: a temp directory to save sam/bam files
# read_dir: contains read files such as :Reads.0.r1.fq.gz
# scripts_code_dir: a directory to save scripts for generating input files 


input_dir="test_data/input" #the dir to the input data
assembly="${input_dir}/final.contigs.fa" #the file of the assembled contigs
mapdir="${input_dir}/map"   # the path to save the alignment files
read_dir="test_data/input/Reads" #the path where you put your reads in
scripts_code_dir="scripts" #the path to our "scripts" files

if [ ! -d $mapdir ]; then
mkdir $mapdir
fi

samtools faidx $assembly;
awk -v OFS='\t' {'print $1,$2'} ${assembly}.fai > ${input_dir}/length.txt;

cnt=0;

for file in ${read_dir}/*;
do echo $file;
let cnt=cnt+1;
echo $cnt;
predix=`basename ${file}`
minimap2 -t 30 -ax sr $assembly $file > ${mapdir}/${predix}.sam;
done


for file in ${mapdir}/*.sam
do
    stub=${file%.sam}
    echo $stub  
    samtools view -h -b -S $file > ${stub}.bam; samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam;samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam; bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g ${input_dir}/length.txt > ${stub}_cov.txt 
done

for i in ${mapdir}/*_cov.txt
do
   echo $i
   stub=${i%_cov.txt}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > ${stub}_cov.csv
done


perl ${scripts_code_dir}/Collate.pl $mapdir > ${input_dir}/coverage.csv

perl -pe "s/,/\t/g;" ${input_dir}/coverage.csv > ${input_dir}/coverage.tsv

