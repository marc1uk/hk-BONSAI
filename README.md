# BONSAI
Low Energy reconstruction code

## How to compile
* Set the environment variable $ROOTSYS to point to your ROOT install
* Set the environment variable $WCSIMDIR to point to your WCSim install
* Set the environment variable $BONSAIDIR to point to this directory
* Run make
  * This will produce libWCSimBonsai.so in $BONSAIDIR

## How to run
* $BONSAIDIR/sample-root-scripts/like.bin should be replaced for each geometry configure.
  * Default like.bin is Cyl60.bin (60m cylinder with 20" B&L PMT.)
  * Other files exist in $BONSAIDIR/sample-root-scripts/likelihood_binary/
    * Cyl100.bin   : 100m cylinder with 20" B&L PMT
    * Cyl80.bin    : 80m cylinder with 20" B&L PMT
    * Cyl60.bin    : 60m cylinder with 20" B&L PMT
    * HK20PMT.bin  : Egg shape HK with 20" Normal PMT
    * HK20BL.bin   : Egg shape HK with 20" B&L PMT
    * SK20BL.bin   : Super-K detector with 20" B&L PMT
    * SK12BL.bin   : Super-K detector with 12" B&L PMT
    * SK20PMT.bin  : Super-K detector with 20" Normal PMT (SK-I ~ IV configure.)
* Run $BONSAIRDIR/sample-root-scripts/sample_bonsai.C on your WCSim output
  * It is recommended to use $BONSAIDIR/rootbonsai (a root wrapper) to run scripts in compiled mode.
    rootbonsai automatically handles telling ROOT of libs and incs
  * Note: currently you need to run in the folder sample_root_scripts (to pickup the .dat and .bin files)

## The important part of the code
* WCSimBonsai is the main class. It reads in geometry information, and calls reconstruction functions.
  It should be run by calling two methods
  * Init(WCSimRootGeom*)
    * Sets up all the variables, and geometry
    * Run once per file
  * BonsaiFit(float *vert,float *result,float *maxlike,int *nsel,
                int *nhit,int *cab,float *t,float *q)
    * float * vert
      * OUTPUT: Reconstructed X,Y,Z,T
    * float * result
      * OUTPUT: ?
    * float * maxlike
      * OUTPUT: [2] is goodness of fit
    * int * nsel
      * ?
    * int * nhit
      * INPUT: Number of hits
    * int * cab
      * INPUT: List of hit tube IDs. Length number of hits
    * float * t
      * INPUT: List of hit times. Length number of hits
    * float * q
      * INPUT: List of hit charges. Length number of hits
    * Run once per trigger
