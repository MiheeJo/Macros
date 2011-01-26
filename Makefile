ROOTCFLAGS	=	$(shell root-config --cflags)
ROOTGLIBS		=	$(shell root-config --glibs)

CPP					=	g++
CPPFLAGS		=	-g -fPIC -Wno-deprecated -O2 -ansi
LD					=	g++
LDFLAGS			=	-g
SOFLAGS			=	-shared

CPPFLAGS		+= $(ROOTCFLAGS)
NGLIBS			=	$(ROOTGLIBS)
#NGLIBS			+= -lMinuit -lRooFit
NGLIBS			+= -L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/cms/cmssw/CMSSW_3_6_3/external/slc5_ia32_gcc434/lib -lMinuit -lRooFit
GLIBS				= $(filter-out -lNew, $(NGLIBS))

CPP         += -I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms6/include
OUTLIB			= ../lib/

.SUFFIXES:	.cc,.C,.hh,.h
.PREFIXES:	../lib/

readDataSet:	./readDataSet.cpp
	$(CPP) $(CPPFLAGS) -o Fit $(OUTLIB)*.o $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o $(OUTLIB)*.so

all:	readDataSet
