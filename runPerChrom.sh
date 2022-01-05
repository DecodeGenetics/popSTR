#!/bin/bash
set -e
set -o pipefail
if [[ "$#" -ne 5 ]]; then
  echo "Usage: runPerChrom.sh <bamList> <reference> <chrom> <markersPerJob> <nCores>"
  exit 1
fi

BAMLIST=$1
REFERENCE=$2
CHROM=$3
MARKERS_PER_JOB=$4
N_CORES=$5
CODE_DIR=`dirname $0`
CURR_DIR=`pwd`

TOTAL_MARKERS=`wc -l ${CODE_DIR}/markerInfo/${CHROM}markerInfo | cut -d ' ' -f 1`
N_JOBS=`expr ${TOTAL_MARKERS} / ${MARKERS_PER_JOB}`
lastJobIdx=`expr ${N_JOBS} + 1`
echo "Chromosome will be split into ${lastJobIdx} runs."

nDone=`expr ${N_JOBS} \* ${MARKERS_PER_JOB}`
tailNum=`expr ${TOTAL_MARKERS} - ${nDone}`

#run computeReadAttributes, each batch of markers at a time
echo "Computing read attributes."
for ((i=1; i<=$lastJobIdx; i++))
do
    if [ $i -eq $lastJobIdx ]
    then
        echo "${CODE_DIR}/popSTR computeReadAttributes ${BAMLIST} ${CURR_DIR} <(tail -n ${tailNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo | cut -d ' ' -f 1-11,14-) 8 135 ${CHROM} ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N"
    else
        headNum=`expr ${i} \* ${MARKERS_PER_JOB}`
        echo "${CODE_DIR}/popSTR computeReadAttributes ${BAMLIST} ${CURR_DIR} <(head -n ${headNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${MARKERS_PER_JOB} | cut -d ' ' -f 1-11,14-) 8 135 ${CHROM} ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N"
    fi
done | parallel -j $N_CORES

#check if pnSlippage has been computed
pnSlippageFile=${CURR_DIR}/pnSlippage
if [ ! -f "$pnSlippageFile" ]; then
    #Check if attributes for kernel markers have been computed
    chr21dir=${CURR_DIR}/attributes/chr21
    if [ ! -d "$chr21dir" ]; then
        echo "Attributes have not been computed for kernel markers and will be computed now."
        ${CODE_DIR}/popSTR computeReadAttributes ${BAMLIST} ${CURR_DIR} <(cut -d ' ' -f 1-11,14- ${CODE_DIR}/kernel/kernelMarkersInfo) 8 135 chr21 ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N
    fi
    echo "pn-slippage rates have not been estimated and will be estimated now."
    ${CODE_DIR}/popSTR computePnSlippageDefault -PL <(awk '{print $1}' $BAMLIST) -AD ${CURR_DIR}/attributes/chr21 -OF pnSlippage -FP 1 -MS ${CODE_DIR}/kernel/kernelSlippageRates -MD ${CODE_DIR}/kernel/kernelModels
fi

#make directory for vcf files
echo "Making directory for vcf file."
mkdir -p vcfs

#run msGenotyper on all chromosomese using attributes and pnSlippage we now have
echo "Genotyping markers."
for ((i=1; i<=$lastJobIdx; i++))
do
    if [ $i -eq $lastJobIdx ]
    then
        echo "${CODE_DIR}/popSTR msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ${CURR_DIR}/vcfs -VN ${CHROM} -ML <(cut -d ' ' -f 1,2,3,4,8,12,13 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${tailNum}) -I ${lastJobIdx} -FP 1"
    else
        headNum=`expr ${i} \* ${MARKERS_PER_JOB}`
        echo "${CODE_DIR}/popSTR msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ${CURR_DIR}/vcfs -VN ${CHROM} -ML <(cut -d ' ' -f 1,2,3,4,8,12,13 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | head -n ${headNum} | tail -n ${MARKERS_PER_JOB}) -I ${i} -FP 1"
    fi
done | parallel -j $N_CORES

#merge vcf files
echo "Merging vcf files"
grep ^# ${CURR_DIR}/vcfs/${CHROM}_1.vcf > ${CURR_DIR}/vcfs/${CHROM}.vcf
nHeaderLines=`wc -l ${CURR_DIR}/vcfs/${CHROM}.vcf | cut -d ' ' -f 1`
((nHeaderLines++))
for ((i=1; i<=$lastJobIdx; i++))
do
    tail -n+${nHeaderLines} ${CURR_DIR}/vcfs/${CHROM}_${i}.vcf >> ${CURR_DIR}/vcfs/${CHROM}.vcf
    rm ${CURR_DIR}/vcfs/${CHROM}_${i}.vcf
done

#gzipping merged bamFile
bgzip ${CURR_DIR}/vcfs/${CHROM}.vcf
tabix ${CURR_DIR}/vcfs/${CHROM}.vcf.gz

#Merging markeSlippageFiles
echo "Merging markerSlippageFiles"
for ((i=1; i<=$lastJobIdx; i++))
do
    cat markerSlippage${CHROM}_${i} >> markerSlippage${CHROM}
    rm markerSlippage${CHROM}_${i}
done
