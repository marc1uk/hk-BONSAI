#ifndef WCSim_Bonsai
#define WCSim_Bonsai

//////////////////////////////////////////////////////////////////////////
//                                                                      
// WCSim_Bonsai                                                      
//                                                                      
// This class contains information needed to be passed to reconstruction
//     routines.  It's just simple right now-- only the bare-bones  
//     WC info
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TClonesArray.h"
#include "WCSimRootGeom.hh"

class TDirectory;

//////////////////////////////////////////////////////////////////////////
class WCSimBonsai : public TObject {

private:

public:

  WCSimBonsai();
  virtual ~WCSimBonsai();

  // Sets and gets

  Int_t Init(WCSimRootGeom *fGeo, bool swapYandZ);
  Int_t BonsaiFit(float *vert,float *result,float *maxlike,int *nsel,
		             int *nhit,int *cab,float *t,float *q);

  ClassDef(WCSimBonsai,1)  //WCSimRootEvent structure
};

#endif
