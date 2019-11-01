# include libraries
PWD:=$(shell pwd)
#BOOST=$(PWD)/boost/1.61.0
BOOST_INCLUDE=/usr/include/boost # Change to your paths
BOOST_LIB=/usr/lib64
HTSLIB=$(PWD)/htslib-develop
SEQAN=$(PWD)/SeqAnHTS

CXXFLAGS+=-I.
CXXFLAGS+=-isystem $(SEQAN)/include
CXXFLAGS+=-I$(HTSLIB)
CXXFLAGS+=-isystem $(BOOST_INCLUDE)
CXXFLAGS+=-pthread 
CXXFLAGS+=-Wfatal-errors

LDFLAGS=-g -L$(HTSLIB) -Wl,-rpath,$(HTSLIB) -lz -lhts -L$(BOOST_LIB) -Wl,-rpath,$(BOOST_LIB) -lboost_iostreams

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1

#Debug build
#CXXFLAGS+=-O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_HAS_ZLIB=1

# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x

all: popSTR
