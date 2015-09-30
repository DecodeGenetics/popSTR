# Use version 4.8.2 of g++
CXX=/opt/rh/devtoolset-2/root/usr/bin/g++

# include SeqAn libraries
-include ../../../libraries/Makefile.inc
CXXFLAGS+=-I/nfs/prog/bioinfo/apps-x86_64/seqan-library/1.4.1/include
CXXFLAGS+=-I/odinn/users/snaedisk



# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -lz -DSEQAN_HAS_ZLIB=1 

# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x

all: msGenotyper computePnSlippage computeReadAttributes getRefSeq
