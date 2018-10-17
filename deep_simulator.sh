#!/bin/bash

if [ $# -lt 2 ]
then
	echo "Usage: ./deep_simulator.sh <input_genome> <sample_read_num> "
	exit
fi



#------- default parameters -----------#
#-> multiprocess
THREAD_NUM=24       #-> this is the thread (or, CPU) number
#-> genome sampling
SAMPLE_MODE=2       #-> choose from the following distribution: 1: beta_distribution, 2: alpha_distribution, 3: mixed_gamma_dis. default: [2]
GENOME_CIRCULAR=0   #-> 0 for NOT circular and 1 for circular. default: [0]
#-> read geneartion
REPEAT_MODE=0.1     #-> value between 0 and 1, 0.1 (default).
#[note]: 0.1 would give the distribution best simulate the real case
#          0 would give distribution whose basecalling result is slightly worse than the real case
#          1 would give the almost perfect basecalling result using Albacore
NOISE_STD=1.0       #-> set the std of random noise of the signal, default as 1.0



#------- input arguments --------------#
FULLFILE=$1
SAMPLE_NUM=$2

#------- init process -----------------#
FILE_N=$(basename "$FULLFILE")
EXTENSION="${FILE_N##*.}"
# use the filename to name the tmp files
FILENAME="${FILE_N%.*}"

NUM=$(fgrep -o '>' $FULLFILE | wc -l)
PREFIX="signal"
PREALI="align"

# the input should be a fasta file
# we should make a tmp directory named after the input file to
# store the tmp files
rm -rf $FILENAME
mkdir -p $FILENAME
cp $FULLFILE $FILENAME/original_genome

# preprocessing, sampling the read
# satisfy the converage and length distritubtion requirement
echo "Executing the preprocessing step..."
python2 sampling_from_genome/sampling.py \
	-i $FILENAME/original_genome \
	-p $FILENAME/sampled_read \
	-n $SAMPLE_NUM \
	-d $SAMPLE_MODE \
	-c $GENOME_CIRCULAR

# pore model translation
# convert the signal to the original range
# signal duplication 
# done within pore model
echo "Finished the preprocessing step!"
echo "Running the context-dependent pore model..."
rm -rf $FILENAME/signal/*
mkdir -p $FILENAME/signal
rm -rf $FILENAME/align/*
mkdir -p $FILENAME/align
source activate tensorflow_cdpm
export DeepSimulatorHome=`pwd`
python2 pore_model/src/main.py \
	-i $FILENAME/sampled_read.fasta \
	-p $FILENAME/signal/$PREFIX \
	-l $FILENAME/align/$PREALI \
	-t $THREAD_NUM  \
	-a $REPEAT_MODE -s $NOISE_STD

# change the signal file to fasta5 file
echo "Finished generate the simulated signals!"
echo "Converting the signal into FAST5 files..."
rm -rf $FILENAME/fast5/*
mkdir -p $FILENAME/fast5
python2 signal_to_fast5/fast5_modify_signal.py \
	-i signal_to_fast5/template.fast5 \
	-s $FILENAME/signal \
	-d $FILENAME/fast5 

# basecalling using albacore
echo "Finished format converting!"
echo "Running Albacore..."
FAST5_DIR="$FILENAME/fast5"
FASTQ_DIR="$FILENAME/fastq"
rm -rf $FASTQ_DIR/*
mkdir -p $FASTQ_DIR
source activate basecall
read_fast5_basecaller.py -i $FAST5_DIR -s $FASTQ_DIR \
	-c r94_450bps_linear.cfg -o fastq -t $THREAD_NUM

# check result
echo "Basecalling finished!"
echo "Checking the read accuracy..."
cat $FILENAME/fastq/workspace/pass/*.fastq > $FILENAME/test.fastq
mapping_check/minimap2 -Hk19 -t $THREAD_NUM -c $FULLFILE \
	$FILENAME/test.fastq > $FILENAME/mapping.paf
accuracy=`awk 'BEGIN{a=0;b=0}{a+=$10/$11;b++}END{print a/b}' $FILENAME/mapping.paf`
echo "Here is the mapping identity: $accuracy"
echo $accuracy > $FILENAME/accuracy

#---------- exit -----------#
exit 0


