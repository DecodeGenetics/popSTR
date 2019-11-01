#!/bin/bash
BAMLIST=$1
REFERENCE=$2
CODE_DIR=`dirname $0`

#Compute attributes for chr21 markers (from kernel)
echo "${CODE_DIR}/popSTR computeReadAttributes ${BAMLIST} . <(cut -d ' ' -f 1-11,14- ${CODE_DIR}/kernel/kernelMarkersInfo) 8 135 chr21 ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N"
${CODE_DIR}/popSTR computeReadAttributes ${BAMLIST} . <(cut -d ' ' -f 1-11,14- ${CODE_DIR}/kernel/kernelMarkersInfo) 8 135 chr21 ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N

#Compute pnSlippage for samples in BAMLIST
echo "${CODE_DIR}/popSTR computePnSlippageDefault -PL <(awk '{print $1}' $BAMLIST) -AD ./attributes/chr21 -OF pnSlippage -FP 1 -MS ${CODE_DIR}/kernel/kernelSlippageRates -MD ${CODE_DIR}/kernel/kernelModels"
${CODE_DIR}/popSTR computePnSlippageDefault -PL <(awk '{print $1}' $BAMLIST) -AD ./attributes/chr21 -OF pnSlippage -FP 1 -MS ${CODE_DIR}/kernel/kernelSlippageRates -MD ${CODE_DIR}/kernel/kernelModels

#make directory for vcf files
echo "mkdir -p vcfs"
mkdir -p vcfs

#run msGenotyper
echo "${CODE_DIR}/popSTR msGenotyperDefault -ADCN ./attributes/chr21 -PNS pnSlippage -MS markerSlippageChr21 -VD ./vcfs -VN chr21_small -ML <(cut -d ' ' -f 1,2,3,4,12,13 ${CODE_DIR}/kernel/kernelMarkersInfo) -I 0 -FP 1"
${CODE_DIR}/popSTR msGenotyperDefault -ADCN ./attributes/chr21 -PNS pnSlippage -MS markerSlippageChr21 -VD ./vcfs -VN chr21_small -ML <(cut -d ' ' -f 1,2,3,4,12,13 ${CODE_DIR}/kernel/kernelMarkersInfo) -I 0 -FP 1
