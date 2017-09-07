# include libraries
CXXFLAGS+=-I.
CXXFLAGS+=-isystem./seqan-library-1.4.1/include
CXXFLAGS+=-isystem /odinn/users/hannese/git/boost_1_61_0/build/include

LDFLAGS=-g -L/odinn/users/hannese/git/boost_1_61_0/build/lib -lz -lboost_iostreams

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1 
LDLIBS=-lz
# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x

all: msGenotyper computePnSlippage computeReadAttributes getRefSeq computePnSlippageDefault msGenotyperDefault
