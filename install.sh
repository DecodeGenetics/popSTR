#!/bin/bash
echo "Unzipping and compiling external libraries"
unzip -qq SeqAnHTS.zip 
unzip -qq htslib.zip
make -C htslib-develop
unzip -qq liblinear-2.01.zip
echo "Compiling binaries"
make
