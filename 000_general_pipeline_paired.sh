#!/bin/bash

#This is a standard pipeline analysis used in the paper: 
#"Do you cov me? Effect of coverage reduction on species identification and 
#genome reconstruction in complex biological matrices by metagenome shotgun high-throughput sequencing" 

#Author: Fabio Marroni

#We assume that you have the following software packages installed (links provided)
# bbduk: https://sourceforge.net/projects/bbmap/
# erne: http://erne.sourceforge.net/
# kraken: https://ccb.jhu.edu/software/kraken2/index.shtml
# (optional) krona and krona tools: https://github.com/marbl/Krona/wiki 
# megahit: https://github.com/voutcn/megahit
# busco: https://busco.ezlab.org/
#We assume that all the software are in your path. If not, you will have to edit this script.


# Reads file (or a symbolic link to them) should be placed in a writeable folder (called "00_reads", but the name is not important). 
# All commands should be executed from the folder carrying the reads. 
# All results will be placed in sub-folders at the same level of the "00_reads" directory. 

#Different formats of the read file name may require slight adjustments of the pipeline.

###################################
#
#  ____  ____  _____  _    _ _  __
# |  _ \|  _ \|  __ \| |  | | |/ /
# | |_) | |_) | |  | | |  | | ' / 
# |  _ <|  _ <| |  | | |  | |  <  
# | |_) | |_) | |__| | |__| | . \ 
# |____/|____/|_____/ \____/|_|\_\
#
###################################

#Step 1 BBDUK
BBMAPDIR=/projects/novabreed/share/software/bbmap
#Paired reads
read1=SRR6287060_1.fastq.gz
read2=SRR6287060_2.fastq.gz
prefix1=${read1/.fastq.gz/}
prefix2=${read2/.fastq.gz/}
trimmed_dir=../trimmed
mkdir -p $trimmed_dir
#Paired-end
${BBMAPDIR}/bbduk.sh in1=${read1} in2=${read2} out1=${trimmed_dir}/${prefix1}_noerne.fastq.gz out2=${trimmed_dir}/${prefix2}_noerne.fastq.gz ref=${BBMAPDIR}/resources/adapters.fa ktrim=r ktrim=l qtrim=rl trimq=10 tpe tbo
#Paired-end in case you also want to remove low-complexity reads
${BBMAPDIR}/bbduk.sh in1=${read1} in2=${read2} out1=${trimmed_dir}/${prefix1}_noerne.fastq.gz out2=${trimmed_dir}/${prefix2}_noerne.fastq.gz outm1=${trimmed_dir}/${prefix1}_discarded_noerne.fastq.gz outm2=${trimmed_dir}/${prefix2}_discarded_noerne.fastq.gz ref=${BBMAPDIR}/resources/adapters.fa ktrim=r ktrim=l qtrim=rl trimq=10 entropy=0.8



###################################
# ______ _____  _   _ ______ 
# |  ____|  __ \| \ | |  ____|
# | |__  | |__) |  \| | |__   
# |  __| |  _  /| . ` |  __|  
# | |____| | \ \| |\  | |____ 
# |______|_|  \_\_| \_|______|
###################################
#Step 2 erne-filter

prefix=${prefix1/_1/}
min_length=50
    erne-filter \
    --query1 ${read1} \
    --query2 ${read2} \
    --min-size ${min_length} \
    --threads 4 \
    --gzip \
    --output-prefix ${trimmed_dir}/${prefix} > ${trimmed_dir}/${prefix}_ernefilt.log

   
###################################
#                  _   _    
#                 | | | |   
#   ___  ___  __ _| |_| | __
#  / __|/ _ \/ _` | __| |/ /
#  \__ \  __/ (_| | |_|   < 
#  |___/\___|\__, |\__|_|\_\
#               | |         
#               |_|            
###################################
#Step 3 seqtk for size reduction (optional)

trimread1=${prefix1}.fastq.gz
trimread2=${prefix2}.fastq.gz
HOWMANY=500000 #How many reads you want to sample
mkdir -p ../sub_${HOWMANY}
REP=1       #Iteration number (to estimate variation of results at the same sampling level you may want to repeat the sampling in several replicate experiments)
seqtk sample -s$REP$HOWMANY $trimmed_dir/${trimread1} $HOWMANY > ../sub_${HOWMANY}/${trimread1/.fastq.gz/.fastq}; gzip -f ../sub_${HOWMANY}/${trimread1/.fastq.gz/.fastq}
seqtk sample -s$REP$HOWMANY $trimmed_dir/${trimread2} $HOWMANY > ../sub_${HOWMANY}/${trimread2/.fastq.gz/.fastq}; gzip -f ../sub_${HOWMANY}/${trimread2/.fastq.gz/.fastq}
              
###################################
# _  _______            _  ________ _   _ 
# | |/ /  __ \     /\   | |/ /  ____| \ | |
# | ' /| |__) |   /  \  | ' /| |__  |  \| |
# |  < |  _  /   / /\ \ |  < |  __| | . ` |
# | . \| | \ \  / ____ \| . \| |____| |\  |
# |_|\_\_|  \_\/_/    \_\_|\_\______|_| \_|
###################################
#Create database
KDB=/projects/igats/metagenomics/custom_databases/kraken_nt_2.0.6
#Step 4.1: Download taxonomy (using interactive nodes)
kraken2-build --download-taxonomy --db $KDB
#Step 4.2: Download and format library
module load aligners/blast/latest; kraken2-build --download-library nt --threads 12 --db $KDB
#Step 4.3: Build the database
module load aligners/blast/latest; kraken2-build --build --threads 12 --db $KDB
#Step 4.4: Clean the database
kraken2-build --clean --db $KDB


#Step 4.5: run kraken
mkdir -p ../sub_${HOWMANY}/kraken2
kraken2 --threads 16 --paired --gzip-compressed --db $KDB ../sub_${HOWMANY}/${prefix1}.fastq.gz ../sub_${HOWMANY}/${prefix2}.fastq.gz --output ../sub_${HOWMANY}/kraken2/${prefix}.kraken --use-names --report ../sub_${HOWMANY}/kraken2/${prefix}.kraken.report.txt

                                          
                                          
###################################
# __  __ ______ _____          _    _ _____ _______ 
#|  \/  |  ____/ ____|   /\   | |  | |_   _|__   __|
#| \  / | |__ | |  __   /  \  | |__| | | |    | |   
#| |\/| |  __|| | |_ | / /\ \ |  __  | | |    | |   
#| |  | | |___| |__| |/ ____ \| |  | |_| |_   | |   
#|_|  |_|______\_____/_/    \_\_|  |_|_____|  |_|   
#
###################################
#Step 5, perform de novo assembly with megahit  
mkdir -p ../sub_${HOWMANY}/megahit
megahit --no-mercy -t 8 -o ../sub_${HOWMANY}/megahit/${prefix} -1 ../trimmed/$trimread1 -2 ../trimmed/$trimread2                  

###################################
# ____  _    _  _____  _____ ____  
#|  _ \| |  | |/ ____|/ ____/ __ \ 
#| |_) | |  | | (___ | |   | |  | |
#|  _ <| |  | |\___ \| |   | |  | |
#| |_) | |__| |____) | |___| |__| |
#|____/ \____/|_____/ \_____\____/ 
#                                   
###################################
#Busco analysis on the reconstructed contigs (metagenome, without genome binning)
BUSCOPATH=/projects/novabreed/share/mvidotto/bioinfotree/binary/appliedgenomics/local/stow/busco-3.02b/
mkdir -p ../sub_${HOWMANY}/busco
cd ../sub_${HOWMANY}/busco
ASSEMBLY=../megahit/${prefix}/final.contigs.fa
run_busco -i $ASSEMBLY -o ${prefix} -c 8 -l $BUSCOPATH/eukaryota_odb9 -m genome 
#cd ../../00_reads
