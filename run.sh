#!/bin/bash
BAMLIST=$1
#run computeReadAttributes for all chromosomes in parallel
for i in {1..22}
do
    echo "computeReadAttributes ${BAMLIST} . markerInfo/chr${i}markerInfo 8 135 chr${i}"
done | parallel

#run computePnSlippageDefault to get pnSlipps using kernel
computePnSlippageDefault -PL <(cols 1 $BAMLIST) -AD ./attributes/chr21 -OF pnSlippage -FP 1 -MS ./kernel/kernelSlippageRates -MD ./kernel/kernelModels

#make directory for vcf files
mkdir vcfs
#run msGenotyper on all chromosomese using attributes and pnSlippage we now have
for i in {1..22}
do
    echo "msGenotyperDefault -ADCN ./attributes/chr${i} -PNS pnSlippage -MS markerSlippageChr${i} -VD ./vcfs -VN chr${i} -ML <(cols 1,2,3,4 markerInfo/chr${i}markerInfo) -I 0 -FP 1"
done | parallel

