#!/bin/bash
echo "Unzipping external libraries"
unzip -qq SeqAnHTS.zip 
unzip -qq htslib-1.9.zip 
unzip -qq boost-1.61.0.zip 
unzip -qq liblinear-2.01.zip
echo "Unzipping and relocating markerInfo files"
mkdir markerInfo
for i in {1..22}; do gunzip -c chr${i}markerInfo.gz > markerInfo/chr${i}markerInfo; rm chr${i}markerInfo.gz; done
echo "Unpacking kernel"
gunzip kernel.tar.gz
tar -xvf kernel.tar
echo "Compiling binaries"
make
