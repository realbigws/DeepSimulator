#!/bin/bash

if [ $# -lt 1 ]
then
	echo "Usage: ./context_poremodel.sh <input_genomic_sequence> "
	exit
fi


#------ read input arguments ---------#
FULLFILE=$1
PREFIX=signal
THREAD_NUM=8
home=`pwd`

#-> get query_name
fulnam=`basename $FULLFILE`
relnam=${fulnam%.*}


#-> call context-dependent pore model
source activate tensorflow_cdpm
export DeepSimulatorHome=$home
python2 $home/pore_model/src/context_poremodel.py \
	-i $FULLFILE \
	-o $PREFIX \
	-t $THREAD_NUM 
	
source deactivate

