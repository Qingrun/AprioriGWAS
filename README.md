# AprioriGWAS
AprioriGWAS is fast tool to identify gene-gene interaction in genome wide association studies

# Requirements
Java Version 8
plink

# Installation
No installation needed. Java is a cross platform,  jar file can be distributed and executed on Windows systems, Mac systems, etc.

# FAQ
What's the input file
How to run AprioriGWAS

AprioriGWAS usage:

Input files can be .csv file or .tped & .tfam format of plink, example data: 

Type in "java -jar AprioriGWAS_2019.jar ", it will show all the steps, 
Type in "java -jar AprioriGWAS_2019.jar " step_name, it will show what parameters are needed 

If input data is .tped format, you will need .tfam file for sample information, please use plink to do data cleaning

Step 1: use plink to do data cleaning

Step 1.a: Run plink to calculate allele frequency, remove individuals with missing data more than 10%, 

example command:
plink --tfile data --mind 0.1 --recode --out data_remove_ind_missing

Step 1.b: check sample relationship: please refer to this page
https://zzz.bwh.harvard.edu/plink/ibdibs.shtml
example command:
plink --file data_remove_ind_missing --genome

follow the instruction remove related individuals



Step 1.c: remove variants with minor allele frequency less than 5%.

plink --file data  --recode --out cleaned


If input data is .csv format, please start from step 3

Step 3: Create HDF5 file to store genotype information

java -jar AprioriGWAS_2020.jar importhdf5 -ig test.csv -o test.hdf5

Step 4: Do single marker test

java -jar AprioriGWAS_2020.jar SingleMarkerTest -g test.hdf5 -o test_single.txt

Step 5: Select SNP for conditional permutation.  The paremeter for -l equals (int)(-log10(The Most Significant pvalue))+1

java -jar AprioriGWAS_2020.jar SNP4Permutaton -s test_single.txt -l 4

Step 6: Prepare permutation

java -jar AprioriGWAS_2020.jar PreparePermute -g test.hdf5  -f 314 -o permutation/3_314 -n 100
java -jar AprioriGWAS_2020.jar PreparePermute -g test.hdf5  -f 322 -o permutation/0_322 -n 100
java -jar AprioriGWAS_2020.jar PreparePermute -g test.hdf5  -f 310 -o permutation/1_310 -n 100
java -jar AprioriGWAS_2020.jar PreparePermute -g test.hdf5  -f 1043 -o permutation/5_1043 -n 1000

Step 7: Do permutation. This step you need to write a for loop to submit jobs
This step, as the significant level of focus SNP single marker test goes up, the p-value threshold (-pt threshold) should get more strict. 

Example perl code:

   #!/usr/bin/perl
   use strict;
   use warnings;
   
   my @a=(0..99);
   for my $i (@a){
           system("java -jar AprioriGWAS_2020.jar Permutation -g test.hdf5 -f 310 -p $i -i permutation/1_310/ -it 3 -pt 0.01 -s 2");
   
   }
   
sbatch   -o 625_0.out -e 625.err   --export snp='625',fd='0_625',pt='0.01' simulate.sh
sbatch   -o 625_1.out -e 625_1.err   --export snp='625',fd='0_625',pt='0.01' simulate1.sh
sbatch   -o 625_2.out -e 625_2.err   --export snp='625',fd='0_625',pt='0.01' simulate2.sh
sbatch   -o 625_3.out -e 625_3.err   --export snp='625',fd='0_625',pt='0.01' simulate3.sh
sbatch   -o 625_4.out -e 625_4.err   --export snp='625',fd='0_625',pt='0.01' simulate4.sh
sbatch   -o 625_5.out -e 625_5.err   --export snp='625',fd='0_625',pt='0.01' simulate5.sh
sbatch   -o 625_6.out -e 625_6.err   --export snp='625',fd='0_625',pt='0.01' simulate6.sh
sbatch   -o 625_7.out -e 625_7.err   --export snp='625',fd='0_625',pt='0.01' simulate7.sh
sbatch   -o 625_8.out -e 625_8.err   --export snp='625',fd='0_625',pt='0.01' simulate8.sh
sbatch   -o 625_9.out -e 625_9.err   --export snp='625',fd='0_625',pt='0.01' simulate9.sh

sbatch   -o 1043_0.out -e 1043_0.err   --export snp='1043',fd='5_1043',pt='0.01' simulate.sh
sbatch   -o 1043_1.out -e 1043_1.err   --export snp='1043',fd='5_1043',pt='0.01' simulate1.sh
sbatch   -o 1043_2.out -e 1043_2.err   --export snp='1043',fd='5_1043',pt='0.01' simulate2.sh
sbatch   -o 1043_3.out -e 1043_3.err   --export snp='1043',fd='5_1043',pt='0.01' simulate3.sh
sbatch   -o 1043_4.out -e 1043_4.err   --export snp='1043',fd='5_1043',pt='0.01' simulate4.sh
sbatch   -o 1043_5.out -e 1043_5.err   --export snp='1043',fd='5_1043',pt='0.01' simulate5.sh
sbatch   -o 1043_6.out -e 1043_6.err   --export snp='1043',fd='5_1043',pt='0.01' simulate6.sh
sbatch   -o 1043_7.out -e 1043_7.err   --export snp='1043',fd='5_1043',pt='0.01' simulate7.sh
sbatch   -o 1043_8.out -e 1043_8.err   --export snp='1043',fd='5_1043',pt='0.01' simulate8.sh
sbatch   -o 1043_9.out -e 1043_9.err   --export snp='1043',fd='5_1043',pt='0.01' simulate9.sh

Step 8: Get threshold

java -jar AprioriGWAS_2020.jar Threshold -f permutation -s 0_625  -p 1000 -l 2
java -jar AprioriGWAS_2020.jar Threshold -f permutation -s 1_974  -p 1000 -l 2
java -jar AprioriGWAS_2020.jar Threshold -f permutation -s 2_1168  -p 1000 -l 2
java -jar AprioriGWAS_2020.jar Threshold -f permutation -s 3_1223 -p 1000 -l 2
java -jar AprioriGWAS_2020.jar Threshold -f permutation -s 4_1203  -p 1000 -l 2
java -jar AprioriGWAS_2020.jar Threshold -f permutation -s 5_1043 -p 1000 -l 2
java -jar AprioriGWAS_2020.jar Threshold -f permutation -s 10_1211 -p 1000 -l 2

Step 9: Get threshold table
since the p-value are sorted from big to small, to choose the top 5% pvalue, one should use 0.95 percentage threshold

java -jar AprioriGWAS_2020.jar Thresholdtable -f permutation -p 1000 -m 10 -n 2 -q 0.05 -o item
java -jar AprioriGWAS_2020.jar Thresholdtable -f permutation -p 1000 -m 10 -n 2 -q 0.95 -o pval

Step 10: Run AprioriGWAS. First create a result folder to store output
Example code as below:

java -jar AprioriGWAS_2020.jar Apriori -g test.hdf5 -l 2 -o result -s 200 -i 0 -pt permutation/pval_0.95_threshold.txt -v test_single.txt 

sbatch -o 0.out -e 0.err --export index='0' apri.sh
sbatch -o 1.out -e 1.err --export index='1' apri.sh
sbatch -o 2.out -e 2.err --export index='2' apri.sh
sbatch -o 3.out -e 3.err --export index='3' apri.sh
sbatch -o 4.out -e 4.err --export index='4' apri.sh
sbatch -o 5.out -e 5.err --export index='5' apri.sh
sbatch -o 6.out -e 6.err --export index='6' apri.sh
	
# Reference
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003627
