#!/bin/bash
BAMLIST=$1
REFERENCE=$2
CHROM=$3
MARKERS_PER_JOB=$4
CODE_DIR=`dirname $0`
CURR_DIR=`pwd`

TOTAL_MARKERS=`wc -l ${CODE_DIR}/markerInfo/${CHROM}markerInfo | cut -d ' ' -f 1`
N_JOBS=`calc ${TOTAL_MARKERS}/${MARKERS_PER_JOB} | cut -d ' ' -f 3 | awk '{printf("%.f\n", $1-0.5)}'`
lastJobIdx=`calc ${N_JOBS}+1 | cut -d ' ' -f 3 | cut -d '.' -f 1`
echo "Chromosome will be split into ${lastJobIdx} runs."

#run computeReadAttributes, each batch of markers at a time
echo "Computing read attributes."
for ((i=1; i<=$N_JOBS; i++))
do
    headNum=`calc ${i}*${MARKERS_PER_JOB} | cut -d ' ' -f 3 | cut -d '.' -f 1`
    echo "${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} <(head -n ${headNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${MARKERS_PER_JOB}) 8 135 ${CHROM} ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N"
done | parallel

nDone=`calc ${N_JOBS}*${MARKERS_PER_JOB} | cut -d ' ' -f 3 | cut -d '.' -f 1`
tailNum=`calc ${TOTAL_MARKERS}-${nDone} | cut -d ' ' -f 3 | cut -d '.' -f 1`
echo "Computing attributes for last ${tailNum} markers"
echo "${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} <(tail -n ${tailNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo) 8 135 ${CHROM} ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N"
${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} <(tail -n ${tailNum} ${CODE_DIR}/markerInfo/${CHROM}markerInfo) 8 135 ${CHROM} ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N

#Check if attributes for kernel markers have been computed
chr21dir=${CURR_DIR}/attributes/chr21
if [ ! -d "$chr21dir" ]; then
    echo "Attributes have not been computed for kernel markers and will be computed now."
    ${CODE_DIR}/computeReadAttributes ${BAMLIST} ${CURR_DIR} ${CODE_DIR}/kernel/kernelMarkersInfo 8 135 chr21 ${REFERENCE} ${CODE_DIR}/markerInfo/longRepeats N
fi

#run computePnSlippageDefault to get pnSlipps using kernel
echo "Computing pn-slippage rates."
${CODE_DIR}/computePnSlippageDefault -PL <(awk '{print $1}' $BAMLIST) -AD ${CURR_DIR}/attributes/chr21 -OF pnSlippage -FP 1 -MS ${CODE_DIR}/kernel/kernelSlippageRates -MD ${CODE_DIR}/kernel/kernelModels

#make directory for vcf files
echo "Making directory for vcf file."
mkdir -p vcfs

#run msGenotyper on all chromosomese using attributes and pnSlippage we now have
echo "Genotyping markers."
for ((i=1; i<=$N_JOBS; i++))
do
    headNum=`calc ${i}*${MARKERS_PER_JOB} | cut -d ' ' -f 3 | cut -d '.' -f 1`
    echo "${CODE_DIR}/msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ${CURR_DIR}/vcfs -VN ${CHROM} -ML <(cut -f 1,2,3,4 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | head -n ${headNum} | tail -n ${MARKERS_PER_JOB}) -I ${i} -FP 1"
done | parallel

echo "Genotyping attributes for last ${tailNum} markers"
echo "${CODE_DIR}/msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ${CURR_DIR}/vcfs -VN ${CHROM} -ML <(cut -f 1,2,3,4 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${tailNum}) -I ${lastJobIdx} -FP 1"
${CODE_DIR}/msGenotyperDefault -ADCN ${CURR_DIR}/attributes/${CHROM} -PNS pnSlippage -MS markerSlippage${CHROM} -VD ${CURR_DIR}/vcfs -VN ${CHROM} -ML <(cut -f 1,2,3,4 ${CODE_DIR}/markerInfo/${CHROM}markerInfo | tail -n ${tailNum}) -I ${lastJobIdx} -FP 1

#merge bam files
echo "Merging BAM files"
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
    cat markerSlippagechr21_${i} >> markerSlippagechr21
    rm markerSlippagechr21_${i}
done
