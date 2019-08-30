#!/bin/bash
BAMLIST=$1
REFERENCE=$2
CODE_DIR=`dirname $0`
CURR_DIR=`pwd`

#run computeReadAttributes for all chromosomes in parallel
echo "Computing read attributes."
for i in {1..22}
do
    echo "${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} ${CODE_DIR}/markerInfo/chr${i}markerInfo 8 135 chr${i} ${REFERENCE}"
done | parallel

#run computePnSlippageDefault to get pnSlipps using kernel
echo "Computing pn-slippage rates."
${CODE_DIR}/computePnSlippageDefault -PL <(cols 1 $BAMLIST) -AD ${CURR_DIR}/attributes/chr21 -OF pnSlippage -FP 1 -MS ${CODE_DIR}/kernel/kernelSlippageRates -MD ${CODE_DIR}/kernel/kernelModels

#make directory for vcf files
echo "Making directory for vcf file."
mkdir -p vcfs

#run msGenotyper on all chromosomese using attributes and pnSlippage we now have
echo "Genotyping markers."
for i in {1..22}
do
    echo "${CODE_DIR}/msGenotyperDefault -ADCN ${CURR_DIR}/attributes/chr${i} -PNS pnSlippage -MS markerSlippageChr${i} -VD ./vcfs -VN chr${i} -ML <(cols 1,2,3,4 ${CODE_DIR}/markerInfo/chr${i}markerInfo) -I 0 -FP 1"
done | parallel
