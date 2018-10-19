#!/bin/bash

# ============== global variables defined here ========= # start
declare FUNC_RET                 #-> for the function echo value
# ============== global variables defined here ========= # end


#---------------------------------------------------------#
##### ===== All functions are defined here ====== #########
#---------------------------------------------------------#


# ----- usage ------ #
function usage()
{
	echo "DeepSimulator v0.11 [Oct-30-2018] "
	echo "    A Deep Learning based Nanopore simulator which can simulate the process of Nanopore sequencing. "
	echo ""
	echo "USAGE:  ./deep_simulator.sh <-i input_genome> [-n simu_read_num] [-o output_root] [-c CPU_num] "
	echo "                    [-C cirular_genome] [-m sample_mode] [-M repeat_mode] [-s noise_std] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i input_genome   : input genome in FASTA format. "
	echo ""
	echo "***** optional arguments *****"
	echo "-n simu_read_num  : the number of reads need to be simulated. [default = 100] "
	echo ""
	echo "-o output_root    : Default output would the current directory. [default = './\${input_name}_DeepSimu'] "
	echo ""
	echo "-c CPU_num        : Number of processors. [default = 8]"
	echo ""
	echo "-C cirular_genome : 0 for linear genome and 1 for circular genome. [default = 0] "
	echo ""
	echo "-m sample_mode    : choose from the following distribution. [default = 2] "
	echo "                    1: beta_distribution, 2: alpha_distribution, 3: mixed_gamma_dis. "
	echo ""
	echo "-M repeat_mode    : a real value between 0 to 1 to determine the length repeat. [default = 0.1] "
	echo "                    '0.1' would give the distribution best simulate the real case, "
	echo "                    '0.0' would give distribution whose basecalling result is slightly worse than the real case, "
	echo "                    '1.0' would give the almost perfect basecalling result using Albacore. "
	echo ""
	echo "-s noise_std      : set the standard deviation (std) of the random noise of the signal. [default = 1.0] "
	echo ""
	echo "-H home           : home directory of DeepSimulator. [default = `pwd`] "
	echo ""
	exit 1
}


#----- check file existence ------# -> this is a function
function file_exist()
{
	#-- input --#
	local file=${1}    #-> 1st input is the file content
	local name=${2}    #-> 2nd input is the file name
	#-- check --#
	if [ -z "$file" ]
	then
		echo "$name is null !!" >&2
		return 1
	fi
	if [ ! -s "$curdir/$file" ]
	then
		if [ ! -s "$file" ]
		then
			echo "$name $file not found !!" >&2
			return 1
		fi
	else
		file=$curdir/$file
	fi
	#-- return --#
	FUNC_RET=$file
	return 0
}

#----- check root existence ------# -> this is a function
function root_exist()
{
	FUNC_RET=""
	#-- input --#
	local root=${1}    #-> 1st input is the root content
	local name=${2}    #-> 2nd input is the root name
	#-- check --#
	if [ -z "$root" ]
	then
		echo "$name is null !!" >&2
		return 1
	fi
	if [ ! -d "$curdir/$root" ]
	then
		if [ ! -d "$root" ]
		then
			echo "$name $root not found !!" >&2
			return 1
		fi
	else
		root=$curdir/$root
	fi
	#-- return --#
	FUNC_RET=$root
	return 0
}


#------------------------------------------------------------#
##### ===== get pwd and check BlastSearchHome ====== #########
#------------------------------------------------------------#

#------ current directory ------#
curdir="$(pwd)"

#-------- check usage -------#
if [ $# -lt 1 ];
then
        usage
fi


#---------------------------------------------------------#
##### ===== All arguments are defined here ====== #########
#---------------------------------------------------------#

#------- required arguments ------------#
FULLFILE=""
output_root=""

#------- optioanl parameters -----------#
SAMPLE_NUM=100      #-> by default, we simulate 100 reads
#-> multiprocess
THREAD_NUM=8        #-> this is the thread (or, CPU) number
#-> genome sampling
SAMPLE_MODE=2       #-> choose from the following distribution: 1: beta_distribution, 2: alpha_distribution, 3: mixed_gamma_dis. default: [2]
GENOME_CIRCULAR=0   #-> 0 for NOT circular and 1 for circular. default: [0]
#-> read geneartion
REPEAT_MODE=0.1     #-> value between 0 and 1, 0.1 (default).
#[note]: 0.1 would give the distribution best simulate the real case
#        0.0 would give distribution whose basecalling result is slightly worse than the real case
#        1.0 would give the almost perfect basecalling result using Albacore
NOISE_STD=1.0       #-> set the std of random noise of the signal, default as 1.0

#------- home directory -----------------#
home=$curdir


#------- parse arguments ---------------#
while getopts ":i:n:o:c:C:m:M:s:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		FULLFILE=$OPTARG
		;;
	#-> optional arguments
	n)
		SAMPLE_NUM=$OPTARG
		;;
	o)
		output_root=$OPTARG
		;;
	c)
		CPU_num=$OPTARG
		;;
	#-> simulator arguments
	C)
		GENOME_CIRCULAR=$OPTARG
		;;
	m)
		SAMPLE_MODE=$OPTARG
		;;
	M)
		REPEAT_MODE=$OPTARG
		;;
	s)
		NOISE_STD=$OPTARG
		;;
	#-> home directory
	H)
		home=$OPTARG
		;;
	#-> default
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done



#---------------------------------------------------------#
##### ===== Part 0: initial argument check ====== #########
#---------------------------------------------------------#

# ------ check home directory ---------- #
if [ ! -d "$home" ]; 
then
	echo "home directory $home not exist " >&2
	exit 1
fi
#-> change to absolute
if [ ! -d "$curdir/$home" ]
then
	if [ ! -d "$home" ]
	then
		echo "home $home not found !!" >&2
		exit 1
	fi
else
	home=$curdir/$home
fi
#echo "home=$home"


#----------- check input genome  -----------#
if [ -z "$FULLFILE" ]
then
	echo "input input_genome is null !!" >&2
	exit 1
fi
file_exist $FULLFILE FULLFILE
retval=$?
if [ "$retval" != "0"  ]
then
	exit 1
else
	FULLFILE=$FUNC_RET
fi
#-> get query_name
fulnam=`basename $FULLFILE`
relnam=${fulnam%.*}


# ------ check output directory -------- #
if [ "$out_root" == "" ]
then
	out_root=$curdir/${relnam}_DeepSimu
fi
dir_out_root=`dirname $out_root`
nam_out_root=`basename $out_root`
if [ "$dir_out_root" == "." ]
then
	if [ "$nam_out_root" == "." ]
	then
		out_root=$curdir
	else
		out_root=$curdir/$nam_out_root
	fi
fi
mkdir -p $out_root
#-> change to absolute
if [ ! -d "$curdir/$out_root" ]
then
	if [ ! -d "$out_root" ]
	then
		echo "outroot $out_root not found !!" >&2
		exit 1
	fi
else
	out_root=$curdir/$out_root
fi


#--------------------------------------------------------#
##### ===== Part 1: DeepSimulator process ====== #########
#--------------------------------------------------------#

#------- init process -----------------#
FILENAME=$out_root
NUM=$(fgrep -o '>' $FULLFILE | wc -l)
PREFIX="signal"
PREALI="align"

# the input should be a fasta file
# we should make a tmp directory named after the input file to
# store the tmp files
cp $FULLFILE $FILENAME/original_genome

# preprocessing, sampling the read
# satisfy the converage and length distritubtion requirement
echo "Executing the preprocessing step..."
python2 $home/sampling_from_genome/sampling.py \
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
export DeepSimulatorHome=$home
python2 $home/pore_model/src/main.py \
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
python2 $home/signal_to_fast5/fast5_modify_signal.py \
	-i $home/signal_to_fast5/template.fast5 \
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
$home/mapping_check/minimap2 -Hk19 -t $THREAD_NUM -c $FULLFILE \
	$FILENAME/test.fastq > $FILENAME/mapping.paf
accuracy=`awk 'BEGIN{a=0;b=0}{a+=$10/$11;b++}END{print a/b}' $FILENAME/mapping.paf`
echo "Here is the mapping identity: $accuracy"
echo $accuracy > $FILENAME/accuracy

#---------- exit -----------#
exit 0



