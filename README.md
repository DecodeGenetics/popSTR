# PopSTR - A Population based microsatellite genotyper

## Downloading and Installing popSTR

1.Download zip of source code: https://github.com/DecodeGenetics/popSTR/archive/master.zip

2.Unzip source code: `unzip popSTR-master.zip`

3.Move to source code directory: `cd popSTR-master/`

4.Unzip library dependencies: `unzip -qq SeqAnHTS.zip` 
                              `unzip -qq htslib-1.9.zip`
                              `unzip -qq boost-1.61.0.zip`
                              `unzip -qq liblinear-2.01.zip`

5.Run make command: `make`

## Running popSTR

PopSTR has three steps and one must iterate between steps 2 and 3 five times !or! (RECOMMENDED) -> use the kernel provided([Description here](#kernelization)): 

### 1. computeReadAttributes - Find useable reads, estimate number of repeats in them and compute their attributes for logistic regression.

Is run for a list of samples over a list of markers.

Call:

    computeReadAttributes bamList outputDirectory markerInfoFile minFlankLength maxRepeatLength chrom

Parameters:

* `bamList` - Two column file listing all samples to be genotyped, first column: A sample identifier. Second column: Path to the bam file. (An index (.bai file) for all bam files must exist in the same directory)
* `outputDirectory` - A folder called attributes will be created in the directory if it doesn't exist. A folder called [chrom] (last parameter above) will be created in the attributes folder if it doesn't exist. One file will be created per marker, located at [outputDirectory]/attributes/[chrom]/[start]_[motif] where start is the start coordinate of the marker and motif is the marker's repeat motif.
* `markerInfoFile` - A file listing the markers to be genotyped, must only contain markers from one chromosome. Marker files for chr1-ch22 are provided. Format for markerInfoFile:

        chrom startCoordinate endCoordinate repeatMotif numOfRepeatsInRef 1000refBasesBeforeStart 1000refBasesAfterEnd repeatSeqFromRef minFlankLeft minFlankRight repeatPurity

* `minFlankLength` - Minimum number of flanking bases required on each side of a repeat for a read to be considered useful.
* `maxRepeatLength` - All alleles with a basePair length above this number will be lumped together into a greater than allele. Should be set close to 0.5 * readLength in IN.bam.
* `chrom` - All markers in markerInfoFile must come from this chromosome.

## Output:
First line contains an offset entry for each PN, to make seeking in the file possible.

The rest of the attribute-file has the following format for each sample specified in the bamList:

1:Sample id 

2: Markerline

    chrom startCoordinate endCoordinate repeatMotif numOfRepeatsInRef numOfReadsFound  A1-initialization A2-initialization

3: Attribute lines (Number of lines is determined by numOfReadsFound from the Markerline):

    numOfRepeatsInRead alignmentQualBefore alignmentQualAfter locationShift mateEditDist repeatPurity ratioOver20inRepeat repeatSeqFromRead 

* numOfRepeatsInRead - The allele reported by the read, tells us how many repeats were found in the read.
* alignmentQualBefore - Quality of realignment on left side of repeat. Ranges from 0 to 1 and higher values indicate more reliable reads.
* alignmentQualAfter - Quality of realignment on right side of repeat. Ranges from 0 to 1 and higher values indicate more reliable reads.
* locationShift - Measures changes from original alignment during the realignment of flanking sequences, common values lie between 0 and 4.
* mateEditDist - Edit distance of aligned base pairs of the mate sequence to the reference, higher values of this attribute are suspicious.
* repeatPurity - Measures how well the repeat sequence from the read matches the repeat sequence from the reference. Ranges from 0 to 1 and higher values indicate more reliable reads.
* ratioOver20inRepeat - Measures the ratio of bases in the repeat sequence with high quality scores. Ranges from 0 to 1 and higher values indicate more reliable reads.
* repeatSeqFromRead - Repeat sequence from the read.

### 2. computePnSlippage - Estimate sample specific slippage rates

Is run for a list of samples over a list of markers.

Call:

    computePnSlippage -ML markerList -PL pnList -AD attributesDirectory -OF outputFile -IN iterationNumber -FP firstPnIdx -NPN minPnsPerMarker -MS markerSlippageFile -MD modelAndLabelDir 
    
Parameters:

* `markerList` - List of markers to estimate the slippage over, must all come from the same chromosome. Format for markerList: 

                    chrom start end motif

* `pnList` - List of samples to estimate slippage for, should only have a single column with same sample identifiers as in column one of bamList used in computeReadAttributes call.
* `attributesDirectory` - Should be the same as the "outputDirectory" parameter in the computeReadAttributes call + /attributes/chrom where chrom is the chromosome of the markers in markerList.
* `outputFile` - Slippage rate will be written to a file with \_iterationNumber appended to this name, i.e. in iteration 0 outputFile\_0 will be created. OBS:if the file already exists, the slippage rates estimated will be appended to it. If iterationNumber > 0 the program will also look for outputFile\_[iterationNumber - 1] to get slippage estimates from previous iteration.
* `iterationNumber` - Zero based number of iteration.
* `firstPnIdx` - Index of the first sample in pnList within the bamList passed to the computeReadAttributes command, to enable faster seeking within the attributes files.
* `minPnsPerMarker` - Minimum number of individuals used to estimate slippage at a given marker so the slippage can be considered reliable. Slippage rates estimated using a small number of samples might not be accurate. 
* `markerSlippageFile` - ONLY applicable if iterationNumber > 0!
                               Should be a path to the markerSlippage file created by msGenotyper in the last iteration.
                               iterationNumber will be appended to the path.
* `modelAndLabelDir` - ONLY applicable if iterationNumber > 0!
                        Should be modelAndLabelDir supplied in the msGenotyper call for this chrom.
                        The program will find within the folder for chrom the label files created in iteration number [iterationNumber]
                        of msGenotyper. It will also locate a file for each marker containing the current logistic regression model for that marker.

### 3. msGenotyper - Determine genotypes, estimate marker slippage rates and train logistic regression models.

Is run for a list of samples(preferably all samples in bamList from computeReadAttributes command) over a list of markers from the same chromosome.

Call:

    msGenotyper -AD attDir -PNS pnList -ML markerList -I intervalIndex -MS markerSlippageFile -MD modelAndLabelDir -IN iterationNumber -FP firstPnIdx -VD vcfOutputDirectory -VN vcfFileName

Parameters:

* `attDir` - A concatenation of the outputDirectory and chrom paratmeters: [outputDirectory]/[chrom] from the computeReadAttributes call of the chromosome being genotyped.
* `pnList` - Should be the outputFile parameter supplied in the computePnSlippage call.
* `markerList` - List of markers to genotype, has only 2 columns: start and motif
* `intervalIndex` - If splitting the chromosome this parameter indicates the index of the interval specified. 
                    When this functionality is not desired the user can specify any number, it will be appended to names of all output files as "\_intervalIndex" and must be manually
                    removed before performing subsequent iterations.
* `markerSlippageFile` - The marker slippage rates estimated will be written to [markerSlippageFile][iterationNumber]\_[intervalIndex] and if iterationNumber>1 the marker slippage rates
                         estimated in the previous iteration will be read from [markerSlippageFile][iterationNumber-1]\_[intervalIndex]
* `modelAndLabelDir` -  In this folder a files called [start]\_[motif][iterationNumber] and model\_[start]\_[motif] will be created for each marker.
                                  If iteration number>1 the labels and models from the previous iteration will be read.
* `iterationNumber` - One based number of iteration.
* `firstPnIdx` - Index of the first sample in pnList within the bamList passed to the computeReadAttributes command, to enable faster seeking within the attributes files.
* `vcfOutputDirectory` - ONLY applicable if iteration number == 5!
                          Directory to place vcf-files containing genotypes for all individuals at all markers contained within the interval being considered.
* `vcfFileName` - ONLY applicable if iteration number == 5! Genotypes will be written to a vcf-file [vcfOutputDirectory]/[vcfFileName]\_[intervalIndex].vcf.
                  
#### Additional notes:

When finishing an iteration of msGenotyper, these steps must be taken before starting the next iteration of computePnSlippage.

1. MarkerSlippage-files for all intervals must be concatenated into a single file for that chromosome.
   For example if one has just finished iteration 1 and the chromosome has been split into 20 intervals, the following command must be run in the directory containing the interval-split
   markerSlippage-files:
   
   `for i in {1..22}; do cat markerSlippageFile1_${i} >> markerSlippageFile1; done`
   
A specific vcf file will be created for each interval and when all intervals have finished running these must be manually concatenated into a single file.

    `grep ^# [vcfOutputDirectory]/vcfFileName_1.vcf > [vcfOutputDirectory]/vcfFileName.vcf`
    `for i in {1..[numIntervals]}; do grep -v ^# [vcfOutputDirectory]/vcfFileName_${i}.vcf >> [vcfOutputDirectory]/vcfFileName.vcf; done`
<a name="kernelization"></a>
#### Kernelization


To circumvent the iterative part using your entire dataset it is possible to train only slippage rates for a selected set of markers (a kernel) and run the "Default" versions of `computePnSlippage` and `msGenotyper` to obtain genotypes directly.
A kernel containing 8779 markers on chr21(HG38 coordinates) is supplied along with the software (kernelSlippageRates, kernelModels(1 file for each marker in the kernel) and kernelMarkersInfo).
To use the kernel one must follow the same three steps as before with slight changes:

##### 1. Run computeReadAttributes as described in the iterative version for all chromosomes on all PNs to be genotyped.

##### 2. Run computePnSlippageDefault for all PNs (the call is different from the one above).

Call:

    `computePnSlippageDefault -PL pnList -AD attributesDirectory -OF outputFile -FP firstPnIdx -MS kernelSlippageRates -MD pathToKernelModelsFiles`

Parameters:

* `pnList` - List of samples to estimate slippage for, should only have a single column with same sample identifiers as in column one of bamList used in computeReadAttributes call.
* `attributesDirectory` - outputDirectory/attributes/chr21 Parameter 2 from running computeReadAttributes with kernelMarkersInfo from the kernel as parameter number 3 appended with attributes and chr21 as shown.
* `outputFile` - Slippage rate will be written(appended if file exists) to this file.
* `firstPnIdx` - Index of the first sample in pnList within the bamList passed to the computeReadAttributes command, to enable faster seeking within the attributes files.
* `kernelSlippageRates` - supplied in the kernel tarball, contains slippage rates and other info for kernel markers.
* `pathToKernelModelsFiles` - Folder supplied in the kernel tarball. Logistic regression models for all markers in the kernel.

##### 3. Run msGenotyperDefault (The call is different from above).

Call:

    `msGenotyperDefault -ADCN attributesDirectory/chromNum/ -PNS pnSlippageFile -MS markerSlippageFile -VD vcfOutputDirectory -VN vcfFileName -ML markerList -I intervalIndex -FP firstPnIdx`

Parameters:

* `attributesDirectory/chromNum/` - Should be a path to output files from running computeReadAttributes in the standard way followed by the name of the chromosome to be genotyped. This directory will be searched for attribute files for all markers in the markerList file.
* `pnSlippageFile` - This file should contain the output of computePnSlippageDefault for all individuals to be genotyped and have 2 columns, sample identifier and slippage rate.
* `markerSlippageFile` - Slippage rates for all markers to be genotyped will be written to a file called [markerSlippageFile]_[intervalIndex].
* `vcfOutputDirectory` - Directory to place vcf-files containing genotypes for all individuals at all markers in the markerList file.
* `vcfFileName` - Genotypes will be written to a vcf-file vcfOutputDirectory/vcfFileName_[intervalIndex].vcf
* `markerList` -  List of markers to estimate the slippage over, must all come from the same chromosome. Format for markerList: 

                    chrom start end motif
                    
* `intervalIndex` - If splitting the chromosome this parameter indicates the index of the interval specified. 
                    When this functionality is not desired the user can specify any number, it will be appended to names of all output files as "_intervalIndex" and must be manually
                    removed before performing subsequent iterations.
* `firstPnIdx` - Index of the first sample in pnSlippageFile within the bamList passed to the computeReadAttributes command, to enable faster seeking within the attributes files.
