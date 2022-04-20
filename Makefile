# include libraries
PWD:=$(shell pwd)
#BOOST=$(PWD)/boost/1.61.0
BOOST_INCLUDE=/usr/include/boost # Change to your paths
BOOST_LIB=/usr/lib64
HTSLIB=$(PWD)/htslib-develop
SEQAN=$(PWD)/SeqAnHTS
OBJS=liblinear.o popSTR.o msGenotyper.o computePnSlippage.o computeReadAttributes.o computePnSlippageDefault.o msGenotyperDefault.o

CXXFLAGS+=-I.
CXXFLAGS+=-isystem $(SEQAN)/include
CXXFLAGS+=-I$(HTSLIB)
CXXFLAGS+=-isystem $(BOOST_INCLUDE)
CXXFLAGS+=-pthread
CXXFLAGS+=-Wfatal-errors

# RELEASE build
CXXFLAGS+=-O3
LDFLAGS=-O3

#Debug build
# CXXFLAGS+=-g -O0
# LDFLAGS=-g -O0

# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_HAS_ZLIB=1 -Wmaybe-uninitialized
LDFLAGS+=-pthread -L$(HTSLIB) -Wl,-rpath,$(HTSLIB) -lz -lhts -L$(BOOST_LIB) -Wl,-rpath,$(BOOST_LIB) -lboost_iostreams

all: popSTR

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

popSTR: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f $(OBJS)
