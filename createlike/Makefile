INCS=-I$(ROOTSYS)/include/ -I.
LIBSROOT=-L$(ROOTSYS)/lib/ -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui
# -lCint

all: showlike create_like

binfile.o: binfile.cc binfile.h
	g++ -c binfile.cc

create_like: create_like.cc binfile.o
	g++ -o create_like -I. create_like.cc binfile.o

showlike: showlike.cc binfile.o
	g++ -std=c++11 -o showlike $(INCS) showlike.cc binfile.o $(LIBSROOT)
