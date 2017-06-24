//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jun 18 09:52:13 2017 by ROOT version 6.06/04
// from TTree cBonsaiEvent/Bonsai Reconstruction Tree
// found on file: bonsaiout.root
//////////////////////////////////////////////////////////

#ifndef cBonsaiEvent_h
#define cBonsaiEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "TVector3.h"
#include "vector"
#include "vector"

class cBonsaiEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   string          wcsimfilestring;
   string*         WCSim_Filep=&wcsimfilestring;
   Int_t           RunId;
   Int_t           EventId;
   Int_t           SubtriggerId;
   Double_t        NumDigits;
   Double_t        TotalTankCharge;
   Double_t        TrueVtx_X;
   Double_t        TrueVtx_Y;
   Double_t        TrueVtx_Z;
   Double_t        TrueVtx_T;
   Double_t        TrueVtx_R;
   TVector3        TrueVtx_Vec;
   TVector3*       TrueVtx_Vecp=&TrueVtx_Vec;
   Int_t           Vertex_Found;
   Bool_t          Vertex_In_Tank;
   Double_t        RecoVtx_X;
   Double_t        RecoVtx_Y;
   Double_t        RecoVtx_Z;
   Double_t        RecoVtx_T;
   Double_t        RecoVtx_R;
   TVector3        RecoVtx_Vec;
   TVector3*       RecoVtx_Vecp=&RecoVtx_Vec;
   Double_t        VtxDiff_X;
   Double_t        VtxDiff_Y;
   Double_t        VtxDiff_Z;
   Double_t        VtxDiff_T;
   TVector3        VtxDiff_Vec;
   TVector3*       VtxDiff_Vecp=&VtxDiff_Vec;
   Double_t        VtxDiff_Perp;
   Double_t        VtxDiff_Par;
   Double_t        VtxDiff_Mag;
   Double_t        Vtx_Goodness;
   Double_t        TrueDir_X;
   Double_t        TrueDir_Y;
   Double_t        TrueDir_Z;
   TVector3        TrueDir_Vec;
   TVector3*       TrueDir_Vecp=&TrueDir_Vec;
   Double_t        RecoDir_X;
   Double_t        RecoDir_Y;
   Double_t        RecoDir_Z;
   Double_t        RecoDir_Zenith;
   Double_t        RecoDir_Azimuth;
   Double_t        RecoDir_Rotation;
   Double_t        RecoDir_CosOpenAng;
   Double_t        RecoDir_Ellipticity;
   TVector3        RecoDir_Vec;
   TVector3*       RecoDir_Vecp=&RecoDir_Vec;
   Double_t        RecoDir_Error;
   Double_t        DirDiff_X;
   Double_t        DirDiff_Y;
   Double_t        DirDiff_Z;
   Double_t        DirDiff_AngMag;
   Double_t        Dir_Goodness;
   Double_t        Fitting_Time;
   vector<double>  NHitsUsedForFit;
   vector<double>* NHitsUsedForFitp=&NHitsUsedForFit;
   TVector3        Tank_Exit_Point;
   TVector3*       Tank_Exit_Pointp=&Tank_Exit_Point;
   Double_t        TrackLengthInTank;
   Double_t        TankEnergyLoss;
   Double_t        TankEnergyLossError;
   Bool_t          Intercepts_MRD;
   TVector3        Projected_MRDEntryPoint;
   TVector3*       Projected_MRDEntryp=&Projected_MRDEntryPoint;
   Int_t           MrdTrackId;

   // List of branches
   TBranch        *b_WCSim_File;   //!
   TBranch        *b_RunId;   //!
   TBranch        *b_EventId;   //!
   TBranch        *b_SubtriggerId;   //!
   TBranch        *b_NumDigits;   //!
   TBranch        *b_TotalTankCharge;   //!
   TBranch        *b_TrueVtx_X;   //!
   TBranch        *b_TrueVtx_Y;   //!
   TBranch        *b_TrueVtx_Z;   //!
   TBranch        *b_TrueVtx_T;   //!
   TBranch        *b_TrueVtx_R;   //!
   TBranch        *b_TrueVtx_Vec;   //!
   TBranch        *b_Vertex_Found;   //!
   TBranch        *b_Vertex_In_Tank;   //!
   TBranch        *b_RecoVtx_X;   //!
   TBranch        *b_RecoVtx_Y;   //!
   TBranch        *b_RecoVtx_Z;   //!
   TBranch        *b_RecoVtx_T;   //!
   TBranch        *b_RecoVtx_R;   //!
   TBranch        *b_RecoVtx_Vec;   //!
   TBranch        *b_VtxDiff_X;   //!
   TBranch        *b_VtxDiff_Y;   //!
   TBranch        *b_VtxDiff_Z;   //!
   TBranch        *b_VtxDiff_T;   //!
   TBranch        *b_VtxDiff_Vec;   //!
   TBranch        *b_VtxDiff_Perp;   //!
   TBranch        *b_VtxDiff_Par;   //!
   TBranch        *b_VtxDiff_Mag;   //!
   TBranch        *b_Vtx_Goodness;   //!
   TBranch        *b_TrueDir_X;   //!
   TBranch        *b_TrueDir_Y;   //!
   TBranch        *b_TrueDir_Z;   //!
   TBranch        *b_TrueDir_Vec;   //!
   TBranch        *b_RecoDir_X;   //!
   TBranch        *b_RecoDir_Y;   //!
   TBranch        *b_RecoDir_Z;   //!
   TBranch        *b_RecoDir_Zenith;   //!
   TBranch        *b_RecoDir_Azimuth;   //!
   TBranch        *b_RecoDir_Rotation;   //!
   TBranch        *b_RecoDir_CosOpenAng;   //!
   TBranch        *b_RecoDir_Ellipticity;   //!
   TBranch        *b_RecoDir_Vec;   //!
   TBranch        *b_RecoDir_Error;   //!
   TBranch        *b_DirDiff_X;   //!
   TBranch        *b_DirDiff_Y;   //!
   TBranch        *b_DirDiff_Z;   //!
   TBranch        *b_DirDiff_AngMag;   //!
   TBranch        *b_Dir_Goodness;   //!
   TBranch        *b_Fitting_Time;   //!
   TBranch        *b_NHitsUsedForFit;   //!
   TBranch        *b_Tank_Exit_Point;   //!
   TBranch        *b_TrackLengthInTank;   //!
   TBranch        *b_TankEnergyLoss;   //!
   TBranch        *b_TankEnergyLossError;   //!
   TBranch        *b_Intercepts_MRD;   //!
   TBranch        *b_Projected_MRDEntry;   //!
   TBranch        *b_MrdTrackId;   //!

   cBonsaiEvent(TTree *tree=0);
   virtual ~cBonsaiEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void             DisableBranches();
};

#endif

#ifdef cBonsaiEvent_cxx
cBonsaiEvent::cBonsaiEvent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("bonsaiout.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("bonsaiout.root");
//      }
//      f->GetObject("bonsaitree",tree);
   } else {
      Init(tree);
   }
}

cBonsaiEvent::~cBonsaiEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t cBonsaiEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t cBonsaiEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void cBonsaiEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

//   // Set object pointer
//   WCSim_Filep = 0;
//   TrueVtx_Vecp = 0;
//   RecoVtx_Vecp = 0;
//   VtxDiff_Vecp = 0;
//   TrueDir_Vecp = 0;
//   RecoDir_Vecp = 0;
//   NHitsUsedForFitp = 0;
//   Tank_Exit_Pointp = 0;
//   Projected_MRDEntryp = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   if(fChain) fChain->ResetBranchAddresses();
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("WCSim_File", &WCSim_Filep, &b_WCSim_File);
   fChain->SetBranchAddress("RunId", &RunId, &b_RunId);
   fChain->SetBranchAddress("EventId", &EventId, &b_EventId);
   fChain->SetBranchAddress("SubtriggerId", &SubtriggerId, &b_SubtriggerId);
   fChain->SetBranchAddress("NumDigits", &NumDigits, &b_NumDigits);
   fChain->SetBranchAddress("TotalTankCharge", &TotalTankCharge, &b_TotalTankCharge);
   fChain->SetBranchAddress("TrueVtx_X", &TrueVtx_X, &b_TrueVtx_X);
   fChain->SetBranchAddress("TrueVtx_Y", &TrueVtx_Y, &b_TrueVtx_Y);
   fChain->SetBranchAddress("TrueVtx_Z", &TrueVtx_Z, &b_TrueVtx_Z);
   fChain->SetBranchAddress("TrueVtx_T", &TrueVtx_T, &b_TrueVtx_T);
   fChain->SetBranchAddress("TrueVtx_R", &TrueVtx_R, &b_TrueVtx_R);
   fChain->SetBranchAddress("TrueVtx_Vec", &TrueVtx_Vecp, &b_TrueVtx_Vec);
   fChain->SetBranchAddress("Vertex_Found", &Vertex_Found, &b_Vertex_Found);
   fChain->SetBranchAddress("Vertex_In_Tank", &Vertex_In_Tank, &b_Vertex_In_Tank);
   fChain->SetBranchAddress("RecoVtx_X", &RecoVtx_X, &b_RecoVtx_X);
   fChain->SetBranchAddress("RecoVtx_Y", &RecoVtx_Y, &b_RecoVtx_Y);
   fChain->SetBranchAddress("RecoVtx_Z", &RecoVtx_Z, &b_RecoVtx_Z);
   fChain->SetBranchAddress("RecoVtx_T", &RecoVtx_T, &b_RecoVtx_T);
   fChain->SetBranchAddress("RecoVtx_R", &RecoVtx_R, &b_RecoVtx_R);
   fChain->SetBranchAddress("RecoVtx_Vec", &RecoVtx_Vecp, &b_RecoVtx_Vec);
   fChain->SetBranchAddress("VtxDiff_X", &VtxDiff_X, &b_VtxDiff_X);
   fChain->SetBranchAddress("VtxDiff_Y", &VtxDiff_Y, &b_VtxDiff_Y);
   fChain->SetBranchAddress("VtxDiff_Z", &VtxDiff_Z, &b_VtxDiff_Z);
   fChain->SetBranchAddress("VtxDiff_T", &VtxDiff_T, &b_VtxDiff_T);
   fChain->SetBranchAddress("VtxDiff_Vec", &VtxDiff_Vecp, &b_VtxDiff_Vec);
   fChain->SetBranchAddress("VtxDiff_Perp", &VtxDiff_Perp, &b_VtxDiff_Perp);
   fChain->SetBranchAddress("VtxDiff_Par", &VtxDiff_Par, &b_VtxDiff_Par);
   fChain->SetBranchAddress("VtxDiff_Mag", &VtxDiff_Mag, &b_VtxDiff_Mag);
   fChain->SetBranchAddress("Vtx_Goodness", &Vtx_Goodness, &b_Vtx_Goodness);
   fChain->SetBranchAddress("TrueDir_X", &TrueDir_X, &b_TrueDir_X);
   fChain->SetBranchAddress("TrueDir_Y", &TrueDir_Y, &b_TrueDir_Y);
   fChain->SetBranchAddress("TrueDir_Z", &TrueDir_Z, &b_TrueDir_Z);
   fChain->SetBranchAddress("TrueDir_Vec", &TrueDir_Vecp, &b_TrueDir_Vec);
   fChain->SetBranchAddress("RecoDir_X", &RecoDir_X, &b_RecoDir_X);
   fChain->SetBranchAddress("RecoDir_Y", &RecoDir_Y, &b_RecoDir_Y);
   fChain->SetBranchAddress("RecoDir_Z", &RecoDir_Z, &b_RecoDir_Z);
   fChain->SetBranchAddress("RecoDir_Zenith", &RecoDir_Zenith, &b_RecoDir_Zenith);
   fChain->SetBranchAddress("RecoDir_Azimuth", &RecoDir_Azimuth, &b_RecoDir_Azimuth);
   fChain->SetBranchAddress("RecoDir_Rotation", &RecoDir_Rotation, &b_RecoDir_Rotation);
   fChain->SetBranchAddress("RecoDir_CosOpenAng", &RecoDir_CosOpenAng, &b_RecoDir_CosOpenAng);
   fChain->SetBranchAddress("RecoDir_Ellipticity", &RecoDir_Ellipticity, &b_RecoDir_Ellipticity);
   fChain->SetBranchAddress("RecoDir_Vec", &RecoDir_Vecp, &b_RecoDir_Vec);
   fChain->SetBranchAddress("RecoDir_Error",&RecoDir_Error,&b_RecoDir_Error);
   fChain->SetBranchAddress("DirDiff_X", &DirDiff_X, &b_DirDiff_X);
   fChain->SetBranchAddress("DirDiff_Y", &DirDiff_Y, &b_DirDiff_Y);
   fChain->SetBranchAddress("DirDiff_Z", &DirDiff_Z, &b_DirDiff_Z);
   fChain->SetBranchAddress("DirDiff_AngMag", &DirDiff_AngMag, &b_DirDiff_AngMag);
   fChain->SetBranchAddress("Dir_Goodness", &Dir_Goodness, &b_Dir_Goodness);
   fChain->SetBranchAddress("Fitting_Time", &Fitting_Time, &b_Fitting_Time);
   fChain->SetBranchAddress("NHitsUsedForFit", &NHitsUsedForFit, &b_NHitsUsedForFit);
   fChain->SetBranchAddress("Tank_Exit_Point", &Tank_Exit_Pointp, &b_Tank_Exit_Point);
   fChain->SetBranchAddress("TrackLengthInTank", &TrackLengthInTank, &b_TrackLengthInTank);
   fChain->SetBranchAddress("TankEnergyLoss", &TankEnergyLoss, &b_TankEnergyLoss);
   fChain->SetBranchAddress("TankEnergyLossError", &TankEnergyLossError, &b_TankEnergyLossError);
   fChain->SetBranchAddress("Intercepts_MRD", &Intercepts_MRD, &b_Intercepts_MRD);
   fChain->SetBranchAddress("Projected_MRDEntry", &Projected_MRDEntryp, &b_Projected_MRDEntry);
   fChain->SetBranchAddress("MrdTrackId", &MrdTrackId, &b_MrdTrackId);
   Notify();
}

void cBonsaiEvent::DisableBranches(){
   // disable majority of branches which aren't used in overall event analysis 
   fChain->SetBranchStatus("*",0); //disable all branches
   // event info so we can match it
   fChain->SetBranchStatus("WCSim_File",1);
   fChain->SetBranchStatus("EventId",1);
   fChain->SetBranchStatus("SubtriggerId",1);
   // reconstructed start vertex info
   fChain->SetBranchStatus("Vertex_Found",1);
   fChain->SetBranchStatus("Vertex_In_Tank",1);
   fChain->SetBranchStatus("RecoVtx_Vec",1);
   fChain->SetBranchStatus("RecoVtx_T",1);
   fChain->SetBranchStatus("RecoDir_Vec",1);
   fChain->SetBranchStatus("RecoDir_Error",1);
   // disable these; they don't seem to be useful.
//   fChain->SetBranchStatus("Vtx_Goodness",1);
//   fChain->SetBranchStatus("Dir_Goodness",1);

   // projected track info for connecting with MRD + veto(?) events
   fChain->SetBranchStatus("Tank_Exit_Point",1);
   fChain->SetBranchStatus("TankEnergyLoss",1);
   fChain->SetBranchStatus("TankEnergyLossError",1);
   fChain->SetBranchStatus("Intercepts_MRD",1);
   fChain->SetBranchStatus("Projected_MRDEntry",1);
   fChain->SetBranchStatus("TrackLengthInTank",1);

   // true info, in case we want it to compare
//   fChain->SetBranchStatus("TrueVtx_Vec",1);
//   fChain->SetBranchStatus("TrueVtx_T",1);
//   fChain->SetBranchStatus("TrueDir_Vec",1);
}

Bool_t cBonsaiEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void cBonsaiEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t cBonsaiEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef cBonsaiEvent_cxx
