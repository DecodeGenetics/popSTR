#PopSTR - A Population based microsatellite genotyper

PopSTR has three steps and one must iterate between steps 2 and 3 five times: 

###1. computeReadAttributes - Find useable reads, estimate number of repeats in them and compute their attributes for logistic regression.

Is run per individual, per chromosome.

Call:

    computeReadAttributes IN.bam outputDirectory markerInfoFile minFlankLength maxRepeatLength PN-id

Parameters:

* `IN.bam` - bam file containing reads from the individual specified in PN-id(the last parameter). An index (.bai file) for the bam file must also exist in the same directory. 
* `outputDirectory` - Directory to write output to, a folder called attributes will be created in the directory if it doesn't exist. A folder with the name of the chromosome being run will be created in the attributes folder if it doesn't exist. The output will be written to a file located at outputDirectory/attributes/chr/PN-id
* `markerInfoFile` - A file listing the markers to be genotyped, must only contain markers from one chromosome. Format for markerInfoFile:

        chrom startCoordinate endCoordinate repeatMotif numOfRepeatsInRef 1000refBasesBeforeStart 1000refBasesAfterEnd repeatSeqFromRef

* `minFlankLength` - Minimum number of flanking bases required on each side of a repeat for a read to be considered useful.
* `maxRepeatLength` - All alleles with a basePair length above this number will be lumped together into a greater than allele.
* `PN-id` - The id of the individual being genotyped.

###2. computePnSlippage - Estimate individual specific slippage rates

Is run per individual.

Call:

    `computePnSlippage attributesDirectory PN-id outputFile iterationNumber minPnsPerMarker [markerSlippageDirectory modelAndLabelDir previousSlippageRate]`
    
Parameters:

* `attributesDirectory` - Should be the same as the "outputDirectory" parameter in the computeReadAttributes call. The program will go through the directory, find the attributes folder, in there it will find the folder for each chromosome and within it the list of useful reads and their attributes for the individual being processed.
* `PN-id` - The id of the individual being genotyped.
* `outputFile` - Slippage rate will be written to a file with the iterationNumber appended to this name, i.e. in iteration 0 outputFile0 will be created. This parameter should be the same for all PNs because the PN-id and slipppage rate will be appended to the file.
* `iterationNumber` - Zero based number of iteration.
* `minPnsPerMarker` - Minimum number of individuals used to estimate slippage at a given marker so the slippage can be considered reliable. Slippage rates estimated using too few individuals might not be accurate. 
* `markerSlippageDirectory` - ONLY applicable if iterationNumber > 0!
                               Should be a path to the marker slippage directory supplied in the msGenotyper call.
                               The program will go through the directory, find the folder for each chromosome and within it the markerSlippage file created in iteration
                               number [iterationNumber] of msGenotyper.
* `modelAndLabelDir` - ONLY applicable if iterationNumber > 0!
                        Should be a path to the logistic regression model and current genotype labels directory supplied in the msGenotyper call.
                        The program will go through the directory, find the folder for each chromosome and within it the labels file created in iteration number [iterationNumber]
                        of msGenotyper. It will also locate a file for each marker containing the current logistic regression model for that marker.
* `previousSlippageRate` - ONLY applicable if iterationNumber > 0!
                        The slippage rate estimated for this PN in iteration number [iterationNumber -1] (can be found by: grep <PN-id> [outputFile][iterationNumber-1] | cut -f 2)

###3. msGenotyper - Determine genotypes, estimate marker slippage rates and train logistic regression models.

Is run per chromosome but for all individuals at once.

Call:

    `msGenotyper attDir PN-slippageFile startCoordinate endCoordinate intervalIndex markerSlippageDir modelAndLabelDir iterationNumber chromName [vcfOutputDirectory vcfFileName]`

Parameters:

* `attDir` - attDir should be the outputDirectory supplied in the computeReadAttributes call and chromNum the chromosome being genotyped.
* `PN-slippageFile` - Should be the outputFile parameter supplied in the computePnSlippage call. Should contain one line per individual to be genotyped.
* `startCoordinate`/`endCoordinate` - If the user wishes to split the genotyping by intervals for paralellization purposes the start and end coordinates should be specified here. 
                      To turn this functionality off specify startCoordinate as 0 and endCoordinate>length(chromosome being genotyped)
* `intervalIndex` - If splitting the chromosome this parameter should indicate the index of the interval specified. 
                    When this functionality is not desired the user can specify any number, it will be appended to names of all output files as "_intervalIndex" and must be manually
                    removed before performing subsequent iterations.
* `markerSlippageDir` - The marker slippage rates estimated will be written to [markerSlippageDir]/[chromName]/markerSlippage[iterationNumber]_[intervalIndex] and if iterationNumber>1 the marker slippage rates
                         estimated in the previous iteration will be read from [markerSlippageDir]/[chromName]/markerSlippageFile[iterationNumber-1]
* `modelAndLabelDir` -  In this directory a file called [PN-id]labels[iterationNumber]_[intervalIndex] will be created for each individual being genotyped and for each marker 
                                  a file called [chromNum]_[starCoordinateOfMarker] containing a logistic regression model for that marker will be created.
                                  If iteration number>1 the labels and models from the previous iteration will be read.
* `iterationNumber` - One based number of iteration.
* `chromName` - Name of chromosome to process.
* `vcfOutputDirectory` - ONLY applicable if iteration number == 5!
                          Directory to place vcf-files containing genotypes for all individuals at all markers contained within the interval being considered.
* `vcfFileName` - ONLY applicable if iteration number == 5!
                  Genotypes will be written to a vcf-file vcfOutputDirectory/vcfFileName_[intervalIndex].vcf
                  
####Additional notes:

When finishing an iteration of msGenotyper, two things must be done before starting the next iteration of computePnSlippage.

1. MarkerSlippage-files for all intervals must be concatenated into a single file for that chromosome.
   For example if one has just finished iteration 1 and the chromosome has been split into 20 intervals, the following command must be run in the directory containing the interval-split
   markerSlippage-files:
   
   `for i in {1..22}; do cat markerSlippageFile1_${i} >> markerSlippageFile1; done`
   
2. All label files for each individual must be concatenated into a single label file.
   Continuing the example from above the following command must be run where we assume the file pn-list contains a list of all PN-ids.
   
   `while read -r PN-id; do for i in {1..20}; do cat ${PN-id}labels1_${i} >> ${PN-id}labels1; done; done`
   
A specific vcf file will be created for each interval and when all intervals have finished running these must be manually concatenated into a single file

    `grep ^# [vcfOutputDirectory]/vcfFileName_1.vcf > [vcfOutputDirectory]/vcfFileName.vcf`
    `for i in {1..[numIntervals]}; do grep -v ^# [vcfOutputDirectory]/vcfFileName_${i}.vcf >> [vcfOutputDirectory]/vcfFileName.vcf; done`
