#!/bin/bash
BAMLIST=$1
REFERENCE=$2
CHROM=$3
MARKERS_PER_JOB=$4
CODE_DIR=`dirname $0`
CURR_DIR=`pwd`

TOTAL_MARKERS=`wc -l ${CODE_DIR}/markerInfo/${CHROM}markerInfo | cols 1`
echo $TOTAL_MARKERS
N_JOBS=`calc ${TOTAL_MARKERS}/${MARKERS_PER_JOB} | cols 3 | awk '{printf("%.f\n", $1-0.5)}'`
echo "Chromosome will be split into ${N_JOBS} runs."

#run computeReadAttributes, each batch of markers at a time
echo "Computing read attributes."
for ((i=1; i<=$N_JOBS; i++))
do
    headNum=`calc ${i}*${MARKERS_PER_JOB} | cols 3 | cut -d '.' -f 1`
    echo "${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} <(head -n ${headNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${MARKERS_PER_JOB}) 8 135 ${CHROM} ${REFERENCE}"
done | parallel

nDone=`calc ${N_JOBS}*${MARKERS_PER_JOB} | cols 3 | cut -d '.' -f 1`
tailNum=`calc ${TOTAL_MARKERS}-${nDone} | cols 3 | cut -d '.' -f 1`
echo "Computing attributes for last ${tailNum} markers"
echo "${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} <(tail -n ${tailNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo) 8 135 ${CHROM} ${REFERENCE}"
${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} <(tail -n ${tailNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo) 8 135 ${CHROM} ${REFERENCE}

#run computePnSlippageDefault to get pnSlipps using kernel
echo "Computing pn-slippage rates."
${CODE_DIR}/computePnSlippageDefault -PL <(cols 1 $BAMLIST) -AD ${CURR_DIR}/attributes/chr21 -OF pnSlippage -FP 1 -MS ${CODE_DIR}/kernel/kernelSlippageRates -MD ${CODE_DIR}/kernel/kernelModels

#make directory for vcf files
echo "Making directory for vcf file."
mkdir -p vcfs

#run msGenotyper on all chromosomese using attributes and pnSlippage we now have
echo "Genotyping markers."
for ((i=1; i<=$N_JOBS; i++))
do
    headNum=`calc ${i}*${MARKERS_PER_JOB} | cols 3 | cut -d '.' -f 1`
    echo "${CODE_DIR}/msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ./vcfs -VN ${CHROM} -ML <(cols 1,2,3,4 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | head -n ${headNum} | tail -n ${MARKERS_PER_JOB}) -I ${i} -FP 1"
done | parallel

echo "Genotyping attributes for last ${tailNum} markers"
echo "${CODE_DIR}/msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ./vcfs -VN ${CHROM} -ML <(cols 1,2,3,4 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${tailNum}) -I ${i} -FP 1"
${CODE_DIR}/msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ./vcfs -VN ${CHROM} -ML <(cols 1,2,3,4 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${tailNum}) -I ${i} -FP 1
