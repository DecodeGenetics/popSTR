#!/bin/bash
BAMLIST=$1

#Compute attributes for chr21 markers (from kernel)
echo "computeReadAttributes ${BAMLIST} . ./kernel/kernelMarkersInfo 8 135 chr21"
computeReadAttributes ${BAMLIST} . ./kernel/kernelMarkersInfo 8 135 chr21

#Compute pnSlippage for samples in BAMLIST
echo "computePnSlippageDefault -PL <(cols 1 ${BAMLIST}) -AD ./attributes/chr21 -OF pnSlippage -FP 1 -MS ./kernel/kernelSlippageRates -MD ./kernel/kernelModels"
computePnSlippageDefault -PL <(cols 1 ${BAMLIST}) -AD ./attributes/chr21 -OF pnSlippage -FP 1 -MS ./kernel/kernelSlippageRates -MD ./kernel/kernelModels

#make directory for vcf files
echo "mkdir -p vcfs"
mkdir -p vcfs

#run msGenotyper
echo "msGenotyperDefault -ADCN ./attributes/chr21 -PNS pnSlippage -MS markerSlippageChr21 -VD ./vcfs -VN chr21_small -ML <(cols 1,2,3,4 markerInfo/chr21markerInfo) -I 0 -FP 1"
msGenotyperDefault -ADCN ./attributes/chr21 -PNS pnSlippage -MS markerSlippageChr21 -VD ./vcfs -VN chr21_small -ML <(cols 1,2,3,4 ./kernel/kernelMarkersInfo) -I 0 -FP 1
