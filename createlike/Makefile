INCS=-I$(ROOTSYS)/include/ -I.
LIBSROOT=-L$(ROOTSYS)/lib/ -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui

all: showlike create_like

binfile.o: binfile.cc binfile.h
	g++ -c binfile.cc

create_like: create_like.cc binfile.o
	g++ -o create_like -I. create_like.cc binfile.o

showlike: showlike.cc binfile.o
	g++ -o showlike $(INCS) showlike.cc binfile.o $(LIBSROOT)
