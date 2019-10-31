#!/bin/bash
echo "Unzipping and compiling external libraries"
unzip -qq SeqAnHTS.zip 
unzip -qq htslib.zip
make -C htslib-develop
unzip -qq liblinear-2.01.zip
echo "Unzipping and relocating markerInfo files"
mkdir markerInfo
for i in {1..22}; do gunzip -c chr${i}markerInfo.gz > markerInfo/chr${i}markerInfo; rm chr${i}markerInfo.gz; done
mv longRepeats markerInfo/
mv defaultModel markerInfo/
echo "Unpacking kernel"
gunzip kernel.tar.gz
tar -xvf kernel.tar
echo "Compiling binaries"
make
