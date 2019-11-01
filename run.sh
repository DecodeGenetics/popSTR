#!/bin/bash
BAMLIST=$1
REFERENCE=$2
MARKERS_PER_JOB=$3
CODE_DIR=`dirname $0`
CURR_DIR=`pwd`

#call runPerChrom for all autosomes sequentially
for i in {1..22}
do
    ${CODE_DIR}/runPerChrom.sh $BAMLIST $REFERENCE chr${i} $MARKERS_PER_JOB
done
