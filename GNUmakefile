# GNUmakefile for libWCSimBonsai 2015/08/12 T.Yano
# --------------------------------------------------------------


ROOTCFLAGS   := $(shell root-config --cflags) -DUSE_ROOT -fPIC
ROOTLIBS     := $(shell root-config --libs)

CPPFLAGS  += -Wno-deprecated 
CPPFLAGS  += -I$(PWD)/bonsai
CPPFLAGS  += -I$(ROOTSYS)/include $(ROOTCFLAGS) 
EXTRALIBS += $(ROOTLIBS)
CXXFLAGS  += -D ConstDirC="\"$(CONSTDIR):$(CONSTDIR)/bonsai:\""
#CXXFLAGS  += ""

CXX = g++

WORKDIR := .
TMPDIR := $(WORKDIR)/tmp/WCSimBonsai


BONSAISO    := libWCSimBonsai.so

#BONSAISRC  := ./src/WCSimBonsaiEvent.cc ./include/WCSimBonsaiEvent.hh ./src/WCSimBonsaiGeom.cc ./include/WCSimBonsaiGeom.hh ./include/WCSimPmtInfo.hh ./include/WCSimBonsaiLinkDef.hh

BONSAISRC := ./bonsai/WCSimBonsai.cc ./bonsai/WCSimRootGeom.cc \
			 ./bonsai/fit_param.cc ./bonsai/pmt_geometry.cc ./bonsai/binfile.cc ./bonsai/bonsai.cc ./bonsai/bonsaifit.cc ./bonsai/plato.cc \
			 ./bonsai/searchgrid.cc ./bonsai/fourhitgrid.cc ./bonsai/centroid.cc ./bonsai/hits.cc ./bonsai/hitsel.cc ./bonsai/goodness.cc \
			 ./bonsai/timefit.cc ./bonsai/likelihood.cc ./bonsai/vertex.cc ./bonsai/pot.cc ./bonsai/tree.cc ./bonsai/vertexfit.cc \
			 ./bonsai/bscalls.cc \
			 ./bonsai/WCSimBonsai.hh \
			 ./bonsai/binfile.h ./bonsai/centroid.h ./bonsai/goodness.h ./bonsai/plato.h ./bonsai/timefit.h \
			 ./bonsai/bonsai.h ./bonsai/fit_param.h ./bonsai/hits.h ./bonsai/pmt_geometry.h ./bonsai/tree.h \
			 ./bonsai/bonsaifit.h ./bonsai/fitquality.h ./bonsai/hitsel.h ./bonsai/pot.h ./bonsai/vertex.h \
			 ./bonsai/bscalls.h ./bonsai/fourhitgrid.h ./bonsai/likelihood.h ./bonsai/searchgrid.h ./bonsai/vertexfit.h \
			 ./bonsai/WCSimRootGeom.hh \
			 ./bonsai/WCSimBonsaiLinkDef.hh

BONSAIOBJS	:= $(WORKDIR)/tmp/WCSimBonsai/fit_param.o $(WORKDIR)/tmp/WCSimBonsai/pmt_geometry.o $(WORKDIR)/tmp/WCSimBonsai/binfile.o $(WORKDIR)/tmp/WCSimBonsai/bonsai.o $(WORKDIR)/tmp/WCSimBonsai/bonsaifit.o $(WORKDIR)/tmp/WCSimBonsai/plato.o \
			   $(WORKDIR)/tmp/WCSimBonsai/searchgrid.o $(WORKDIR)/tmp/WCSimBonsai/fourhitgrid.o $(WORKDIR)/tmp/WCSimBonsai/centroid.o $(WORKDIR)/tmp/WCSimBonsai/hits.o $(WORKDIR)/tmp/WCSimBonsai/hitsel.o $(WORKDIR)/tmp/WCSimBonsai/goodness.o \
			   $(WORKDIR)/tmp/WCSimBonsai/timefit.o $(WORKDIR)/tmp/WCSimBonsai/likelihood.o $(WORKDIR)/tmp/WCSimBonsai/vertex.o $(WORKDIR)/tmp/WCSimBonsai/pot.o $(WORKDIR)/tmp/WCSimBonsai/tree.o $(WORKDIR)/tmp/WCSimBonsai/vertexfit.o \
			   $(WORKDIR)/tmp/WCSimBonsai/bscalls.o \
			   $(WORKDIR)/tmp/WCSimBonsai/WCSimBonsai.o \
			   $(WORKDIR)/tmp/WCSimBonsai/WCSimRootGeom.o \
			   $(WORKDIR)/tmp/WCSimBonsai/WCSimBonsaiDict.o





.PHONY: directories

all: directories libWCSimBonsai.so

directories: $(TMPDIR)

$(TMPDIR) :
	mkdir -p $(TMPDIR)


libWCSimBonsai.so : $(BONSAIOBJS) 
	@if [ ! -d $(TMPDIR) ] ; then mkdir $(TMPDIR) ; echo mkdir $(TMPDIR) ;fi
	@$(CXX) -shared -O $^ -o $(BONSAISO) $(BONSAILIBS) $(EXTRALIBS)

./bonsai/WCSimBonsaiDict.cc : $(BONSAISRC)
	@echo Compiling rootcint ...
	rootcint  -f ./bonsai/WCSimBonsaiDict.cc -c -I./bonsai -I$(shell root-config --incdir) \
		WCSimBonsai.hh \
		WCSimRootGeom.hh \
		WCSimBonsaiLinkDef.hh
#		binfile.h centroid.h goodness.h plato.h timefit.h \
#		bonsai.h fit_param.h hits.h pmt_geometry.h tree.h \
#		bonsaifit.h fitquality.h hitsel.h pot.h vertex.h \
#		bscalls.h fourhitgrid.h likelihood.h searchgrid.h vertexfit.h \

rootcint: ./bonsai/WCSimBonsaiDict.cc

$(TMPDIR)/%.o: bonsai/%.cc
	@echo Compiling $*.cc ...
	@if [ ! -d $(TMPDIR) ] ; then mkdir $(TMPDIR) ; echo mkdir $(TMPDIR) ;fi
	@echo $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(TMPDIR)/$(*F).o $<
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $(TMPDIR)/$(*F).o $<

clean :
	@rm -f $(TMPDIR)/*.o
	@rm -f libWCSimBonsai.so
	@rm -f ./bonsai/WCSimBonsaiDict.cc
