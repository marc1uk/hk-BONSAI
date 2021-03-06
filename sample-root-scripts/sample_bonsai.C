//C++
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//ROOT
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TH1.h>
#include <TRint.h>
#include <TColor.h>
#include <TString.h>
#include <math.h>
//#include "../gitver/Bonsai_v0/bonsai/pmt_geometry.h"
//#include "../gitver/Bonsai_v0/bonsai/likelihood.h"
//#include "../gitver/Bonsai_v0/bonsai/goodness.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
//WCSim
#include "../../../WCSim/gitver/wcsim/include/WCSimRootEvent.hh"
#include "../../../WCSim/gitver/wcsim/include/WCSimRootGeom.hh"
//BONSAI
#include "WCSimBonsai.hh"
#endif

#ifndef VERBOSE
//#define VERBOSE 1
#endif

// Simple example of reading a generated Root file
int sample_bonsai(const char *filename=NULL, bool verbose=false)
{
#if !defined(__MAKECINT__)
  // Load the library with class dictionary info
  // (create with "gmake shared")
  char* wcsimdirenv;
  wcsimdirenv = getenv ("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/wcsim/libWCSimRoot.so");
  }else{
    gSystem->Load("../libWCSimRoot.so");
  }
  char* bonsaidirenv;
  bonsaidirenv = getenv ("BONSAIDIR");
  if(bonsaidirenv !=  NULL){
    gSystem->Load("${BONSAIDIR}/libWCSimBonsai.so");
  }else{
    gSystem->Load("../libWCSimBonsai.so");
  }
#endif
  
//  ===============================================================
//  ANNIE WCSim variables: tank -152<X<152, -212<Y<183, 15<Z<320 cm
//  ===============================================================
  int track_pdg=11;   // ******* SET THIS: NEED TO TELL BONSAI TO SEARCH FOR CORRECT PARTICLE TYPE, mu-, e- ...
  bool DrawPlots=false;
  bool isANNIE=true;
  const Float_t tank_start = 15.70;          // front face of the tank in cm
  const Float_t tank_radius = 152.4;         // tank radius in cm
  const Float_t tank_halfheight = 198.;      // tank half height in cm
  const Float_t tank_yoffset = -14.46;       // tank y offset in cm
  
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
  TFile* fileout = new TFile("bonsaiout.root","RECREATE");
  TTree* treeout = new TTree("bonsaitreeout","Bonsai Reconstruction Tree");
  double numdigitsthisevent; // TODO; give this a proper name in the tree/variable. ncherenkovdigihits.
  TBranch* bNumDigits = treeout->Branch("NumDigits",&numdigitsthisevent);
  double totaltankcharge;
  TBranch* bSumQ = treeout->Branch("TotalTankCharge",&totaltankcharge);
  double True_Vertex_X, True_Vertex_Y, True_Vertex_Z, True_Vertex_T, True_Vertex_R;
  TBranch* bTrueVtx_X = treeout->Branch("TrueVtx_X",&True_Vertex_X);
  TBranch* bTrueVtx_Y = treeout->Branch("TrueVtx_Y",&True_Vertex_Y);
  TBranch* bTrueVtx_Z = treeout->Branch("TrueVtx_Z",&True_Vertex_Z);
  TBranch* bTrueVtx_R = treeout->Branch("TrueVtx_T",&True_Vertex_T);
  TBranch* bTrueVtx_T = treeout->Branch("TrueVtx_R",&True_Vertex_R);
  TVector3 True_Vertex_Vec;
  TBranch* bTrueVtx_Vec = treeout->Branch("TrueVtx_Vec",&True_Vertex_Vec);
  double Reco_Vertex_X, Reco_Vertex_Y, Reco_Vertex_Z, Reco_Vertex_T, Reco_Vertex_R;
  TBranch* bRecoVtx_X = treeout->Branch("RecoVtx_X",&Reco_Vertex_X);
  TBranch* bRecoVtx_Y = treeout->Branch("RecoVtx_Y",&Reco_Vertex_Y);
  TBranch* bRecoVtx_Z = treeout->Branch("RecoVtx_Z",&Reco_Vertex_Z);
  TBranch* bRecoVtx_T = treeout->Branch("RecoVtx_Z",&Reco_Vertex_T);
  TBranch* bRecoVtx_R = treeout->Branch("RecoVtx_R",&Reco_Vertex_R);
  TVector3 Reco_Vertex_Vec;
  TBranch* bRecoVtx_Vec = treeout->Branch("RecoVtx_Vec",&Reco_Vertex_Vec);
  double Vertex_Diff_X, Vertex_Diff_Y, Vertex_Diff_Z, Vertex_Diff_T, Vertex_Diff_Perp, Vertex_Diff_Par, Vertex_Diff_Mag;
  TBranch* bVtxDiff_X = treeout->Branch("VtxDiff_X",&Vertex_Diff_X);
  TBranch* bVtxDiff_Y = treeout->Branch("VtxDiff_Y",&Vertex_Diff_Y);
  TBranch* bVtxDiff_Z = treeout->Branch("VtxDiff_Z",&Vertex_Diff_Z);
  TBranch* bVtxDiff_T = treeout->Branch("VtxDiff_T",&Vertex_Diff_T);
  TVector3 Vertex_Diff_Vec;
  TBranch* bVtxDiff_Vec = treeout->Branch("VtxDiff_Vec",&Vertex_Diff_Vec);
  TBranch* bVtxDiff_Perp = treeout->Branch("VtxDiff_Perp",&Vertex_Diff_Perp);
  TBranch* bVtxDiff_Par = treeout->Branch("VtxDiff_Par",&Vertex_Diff_Par);
  TBranch* bVtxDiff_Mag = treeout->Branch("VtxDiff_Mag",&Vertex_Diff_Mag);
  double True_Dir_X, True_Dir_Y, True_Dir_Z;
  TBranch* bTrueDir_X = treeout->Branch("TrueDir_X",&True_Dir_X);
  TBranch* bTrueDir_Y = treeout->Branch("TrueDir_Y",&True_Dir_Y);
  TBranch* bTrueDir_Z = treeout->Branch("TrueDir_Z",&True_Dir_Z);
  TVector3 True_Dir_Vec;
  TBranch* bTrueDir_Vec = treeout->Branch("TrueDir_Vec",&True_Dir_Vec);
  double Reco_Dir_X, Reco_Dir_Y, Reco_Dir_Z, 
    Reco_Dir_Zenith, Reco_Dir_Azimuth, Reco_Dir_Axial, Reco_Dir_CosOpAng, Reco_Dir_Ellip;
  TBranch* bRecoDir_X = treeout->Branch("RecoDir_X",&Reco_Dir_X);
  TBranch* bRecoDir_Y = treeout->Branch("RecoDir_Y",&Reco_Dir_Y);
  TBranch* bRecoDir_Z = treeout->Branch("RecoDir_Z",&Reco_Dir_Z);
  TBranch* bRecoDir_Zenith = treeout->Branch("RecoDir_Zenith",&Reco_Dir_Zenith);
  TBranch* bRecoDir_Azimuth = treeout->Branch("RecoDir_Azimuth",&Reco_Dir_Azimuth);
  TBranch* bRecoDir_Rotation = treeout->Branch("RecoDir_Rotation",&Reco_Dir_Axial);
  TBranch* bRecoDir_CosOpenAng = treeout->Branch("RecoDir_CosOpenAng",&Reco_Dir_CosOpAng);
  TBranch* bRecoDir_Ellipticity = treeout->Branch("RecoDir_Ellipticity",&Reco_Dir_Ellip);
  TVector3 Reco_Dir_Vec;
  TBranch* bRecoDir_Vec = treeout->Branch("RecoDir_Vec",&Reco_Dir_Vec);
  double Dir_Diff_X, Dir_Diff_Y, Dir_Diff_Z, Dir_Diff_AngMag;
  TBranch* bDirDiff_X = treeout->Branch("DirDiff_X",&Dir_Diff_X);
  TBranch* bDirDiff_Y = treeout->Branch("DirDiff_Y",&Dir_Diff_Y);
  TBranch* bDirDiff_Z = treeout->Branch("DirDiff_Z",&Dir_Diff_Z);
  TBranch* bDirDiff_AngMag = treeout->Branch("DirDiff_AngMag",&Dir_Diff_AngMag);
  double vtxgoodness, dirgoodness;
  TBranch* bVtx_Goodness = treeout->Branch("Vtx_Goodness",&vtxgoodness);
  TBranch* bDir_Goodness = treeout->Branch("Dir_Goodness",&dirgoodness);
  
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
  // Position reconstruction
  Float_t histoextent;
  (isANNIE) ? histoextent=400 : histoextent = 2000;
  TH1F *hTrueVtx_X = new TH1F("True VTX_X", "True VTX_X", 100, -histoextent, histoextent);
  TH1F *hTrueVtx_Y = new TH1F("True VTX_Y", "True VTX_Y", 100, -histoextent, histoextent);
  TH1F *hTrueVtx_Z = new TH1F("True VTX_Z", "True VTX_Z", 100, -histoextent, histoextent);
  TH1F *hTrueVtx_T = new TH1F("True VTX_T", "True VTX_T", 200, -1500, 1500);
  TH1F *hTrueVtx_R = new TH1F("True VTX_R", "True VTX_R", 100, 0, histoextent);
  TH1F *hRecoVtx_X = new TH1F("Reconstructed VTX_X", "Reconstructed VTX_X", 100, -histoextent, histoextent);
  TH1F *hRecoVtx_Y = new TH1F("Reconstructed VTX_Y", "Reconstructed VTX_Y", 100, -histoextent, histoextent);
  TH1F *hRecoVtx_Z = new TH1F("Reconstructed VTX_Z", "Reconstructed VTX_Z", 100, -histoextent, histoextent);
  TH1F *hRecoVtx_T = new TH1F("Reconstructed VTX_T", "Reconstructed VTX_T", 100, -100, 100);
  TH1F *hRecoVtx_R = new TH1F("Reconstructed VTX_R", "Reconstructed VTX_R", 100, 0, histoextent);
  TH1F *hVertexGoodness = new TH1F("Bonsai Vertex Fit Goodness", "Bonsai Vertex Fit Goodness", 100, 0, 1);
  TH1F *hVtxDiff_X = new TH1F("Vertex Error X", "Vertex Error X", 100, -200, 200);
  TH1F *hVtxDiff_Y = new TH1F("Vertex Error Y", "Vertex Error Y", 100, -200, 200);
  TH1F *hVtxDiff_Z = new TH1F("Vertex Error Z", "Vertex Error Z", 100, -200, 200);
  TH1F *hVtxDiff_T = new TH1F("Vertex Error T", "Vertex Error T", 100, -200, 200);
  TH1F *hVtxDiff_Tot = new TH1F("Total Vertex Error", "Total Vertex Error", 100, 0, 300);
  TH1F *hVtxDiff_Par = new TH1F("Vertex Error in Parallel Dir","Vertex Error in Parallel Dir", 100, 0, 300);
  TH1F *hVtxDiff_Perp = new TH1F("Vertex Error in Perpendicular Dir","Vertex Error in Perpendicular Dir", 100, 0, 300);
  TH1F *hVtxDiff_Intg = new TH1F("Integrated Total Vertex Error","Integrated Total Vertex Error", 100, 0, 300);
  
  // Direction reconstruction
  TH1F *hTrueDir_X = new TH1F("True X Dir", "True X Dir", 100, -1.5, 1.5);
  TH1F *hTrueDir_Y = new TH1F("True Y Dir", "True Y Dir", 100, -1.5, 1.5);
  TH1F *hTrueDir_Z = new TH1F("True Z Dir", "True Z Dir", 100, -1.5, 1.5);
  TH1F *hRecoDir_X = new TH1F("Reconstructed X Dir", "Reconstructed X Dir", 100, -TMath::Pi(), TMath::Pi());
  TH1F *hRecoDir_Y = new TH1F("Reconstructed Y Dir", "Reconstructed Y Dir", 100, -TMath::Pi(), TMath::Pi());
  TH1F *hRecoDir_Z = new TH1F("Reconstructed Z Dir", "Reconstructed Z Dir", 100, -TMath::Pi(), TMath::Pi());
  TH1F *hDirDiff_X = new TH1F("Reconstructed Direction Error X", "Reconstructed Direction Error X", 100, -1.5, 1.5);
  TH1F *hDirDiff_Y = new TH1F("Reconstructed Direction Error Y", "Reconstructed Direction Error Y", 100, -1.5, 1.5);
  TH1F *hDirDiff_Z = new TH1F("Reconstructed Direction Error Z", "Reconstructed Direction Error Z", 100, -1.5, 1.5);
  TH1F *hDirDiff_AngMag = new TH1F("Reconstructed Direction Error Tot", "Reconstructed Direction Error Tot", 100, 0., 3.);
  TH1F *hDirDiff_Intg = new TH1F("Integrated Total Direction Error", "Integrated Total Direction Error", 100, 0., 3.);
  TH1F *hDirectionGoodness = new TH1F("Bonsai Direction Fit Goodness", "Bonsai Direction Fit Goodness", 100, -10, 10);
  
  Float_t maxdigits, maxcharge;
  (isANNIE) ? maxdigits=200 : maxdigits=2000;
  (isANNIE) ? maxcharge=50 : maxcharge=500;
  TH1F *hNumDigits = new TH1F("Num PMT Digits", "Num PMT Digits", 100, 0, maxdigits);
  TH1F *hNumHitsPerDigit = new TH1F("Num Hits Per Digit", "Num Hits Per Digit", 100, 0, maxdigits);
  TH1F *hNumNoiseOnlyDigits = new TH1F("Num Noise Only Digits", "Num Noise Only Digits", 100, 0, maxdigits);
  TH1F *hNumPartialNoiseDigits = new TH1F("Num Partial Noise Digits", "Num Partial Noise Digits", 100, 0, maxdigits);
  TH1F *hNoiseCharge = new TH1F("Charge from Noise Hits", "Charge from Noise Hits", 100,0,maxcharge);
  TH1F *hPhotonCharge = new TH1F("Charge from Photon Hits", "Charge from Photon Hits", 100,0,maxcharge);
  TH1F *hTotalCharge = new  TH1F("Total Digit Charge", "Total Digit Charge", 100,0,maxcharge);
  //TH1F *hRecoVtx_Y = new TH1F("SumQ", "SumQ", 200, 0, 200);
  
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
  
  gROOT->cd();
  
  WCSimBonsai* bonsai = new WCSimBonsai();
  
  TFile *file;
  // Open the file
  if (filename==NULL){
    file = new TFile("../wcsim.root","read");
  }else{
    file = new TFile(filename,"read");
  }
  if (!file->IsOpen()){
    std::cout << "Error, could not open input file: " << filename << std::endl;
    return -1;
  }
  
  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");
  
  // Get the number of events
  int nevent = tree->GetEntries();
  if(verbose) printf("nevent %d\n",nevent);
  
  // Create a WCSimRootEvent to put stuff from the tree in
  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
  
  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  
  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
  
  
  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  WCSimRootGeom *geo = 0; 
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
    exit(9);
  }
  geotree->GetEntry(0);
  bonsai->Init(geo);
  int numpmts = geo->GetWCNumPMT();
  
  // start with the main "subevent", as it contains most of the info and always exists.
  WCSimRootTrigger* wcsimrootevent;
  
  // Now loop over events
  //nevent=100;
  int numsuccessfulfits=0, numfailedfits=0, numskippedevents=0;
#ifdef VERBOSE
  cout<<nevent<<" events to analyze"<<endl;
#endif
  for (int ev=0; ev<nevent; ev++){
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(verbose){
      printf("********************************************************");
      printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
                wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n", wcsimrootsuperevent->GetNumberOfSubEvents());
      printf("Jmu %d\n", wcsimrootevent->GetJmu());
      printf("Npar %d\n", wcsimrootevent->GetNpar());
      printf("Ntrack %d\n", wcsimrootevent->GetNtrack());
    }
    
    // load true vertex and initial momentum direction from WCSimRootTrack. Need to find the muon.
    Double_t Emax=0; //we'll pick the highest energy muon, we can only compare one vertex to the fit result
    Int_t mutrackindex=-1;
    int numtracks= wcsimrootevent->GetNtrack();
    WCSimRootTrack* nextrack;
#ifdef VERBOSE
    cout<<"find mu"<<endl;
#endif
    // scan through the truth tracks, find the primary muon and pull info from it
    for(int track=0; track<numtracks; track++){
      nextrack = (WCSimRootTrack*)wcsimrootevent->GetTracks()->At(track);
//          a WCSimRootTrack has methods: 
//          Int_t     GetIpnu()             pdg
//          Int_t     GetFlag()             -1: neutrino primary, -2: neutrino target, 0: other
//          Float_t   GetM()                mass
//          Float_t   GetP()                momentum magnitude
//          Float_t   GetE()                energy (inc rest mass^2)
//          Int_t     GetStartvol()         starting volume: 10 is tank, 20 is facc, 30 is mrd
//          Int_t     GetStopvol()          stopping volume: but these may not be set.
//          Float_t   GetDir(Int_t i=0)     momentum unit vector
//          Float_t   GetPdir(Int_t i=0)    momentum vector
//          Float_t   GetStop(Int_t i=0)    stopping vertex x,y,z for i=0-2, in cm
//          Float_t   GetStart(Int_t i=0)   starting vertex x,y,z for i=0-2, in cm
//          Int_t     GetParenttype()       parent pdg, 0 for primary.
//          Float_t   GetTime()             starting time of particle
//          Int_t     GetId()               wcsim trackid
      if(/*nextrack->GetParenttype()==0&&*/nextrack->GetIpnu()==track_pdg&&nextrack->GetE()>Emax){
        Emax=nextrack->GetE();
        mutrackindex=track;
      }
    }
#ifdef VERBOSE
    cout<<"mu found"<<endl;
#endif
    if(mutrackindex<0) {continue; numskippedevents++;}// this event had no muons!
    cout<<"✩ ";
    if(verbose){
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0), 
         wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    } else {
      //std::cout<<"Event "<<ev<<std::endl;
    }
    
    nextrack = (WCSimRootTrack*)wcsimrootevent->GetTracks()->At(mutrackindex);
    if(isANNIE){
      True_Vertex_X=nextrack->GetStart(0);
      True_Vertex_Y=(nextrack->GetStart(2))-tank_start-tank_radius;
      True_Vertex_Z=(nextrack->GetStart(1))-tank_yoffset;
      True_Vertex_R=sqrt(pow(True_Vertex_X, 2) + pow(True_Vertex_Z, 2));
      True_Vertex_T=nextrack->GetTime();
      True_Vertex_Vec=TVector3(True_Vertex_X,True_Vertex_Y,True_Vertex_Z);
      
      True_Dir_X=nextrack->GetDir(0);
      True_Dir_Y=nextrack->GetDir(2);
      True_Dir_Z=nextrack->GetDir(1);
      True_Dir_Vec = TVector3(True_Dir_X,True_Dir_Y,True_Dir_Z);
      
    } else {
    
      True_Vertex_X=nextrack->GetStart(0);
      True_Vertex_Y=nextrack->GetStart(1);
      True_Vertex_Z=nextrack->GetStart(2);
      True_Vertex_R=sqrt(pow(True_Vertex_X, 2) + pow(True_Vertex_Z, 2));
      True_Vertex_T=nextrack->GetTime();
      True_Vertex_Vec=TVector3(True_Vertex_X,True_Vertex_Y,True_Vertex_Z);
      
      True_Dir_X=nextrack->GetDir(0);
      True_Dir_Y=nextrack->GetDir(1);
      True_Dir_Z=nextrack->GetDir(2);
      True_Dir_Vec = TVector3(True_Dir_X,True_Dir_Y,True_Dir_Z);
    }
//    don't fill histograms with true vertex until we've checked if we managed a fit
    
    // Loop over sub events and perform vertex fit
    // ===========================================
    float bsT[500],bsQ[500];
    float bsvertex[4],bsresult[6];
    float bsgood[500];
    int bsCAB[500];
    int bsnhit[1]; //nsel (SLE)
    int bsnsel[2]; //nsel (SLE)
    if(verbose) std::cout << "DIGITIZED HITS:" << std::endl;
    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
    for (int index = 0 ; index < wcsimrootsuperevent->GetNumberOfEvents(); index++){
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
      if(verbose) std::cout << "Sub event number = " << index << "\n";
      
      int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
      if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      bsnhit[0] = ncherenkovdigihits;
#ifdef VERBOSE
      cout<<"filling array of Ts and Qs"<<endl;
#endif
      for (int i=0;i<ncherenkovdigihits;i++){
          // Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
          TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
          
          WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
            dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
          
          bsT[i]=wcsimrootcherenkovdigihit->GetT();
          bsQ[i]=wcsimrootcherenkovdigihit->GetQ();
          bsCAB[i]=wcsimrootcherenkovdigihit->GetTubeId();
          
          if(verbose){
            if ( i < 10 ) // Only print first XX=10 tubes
              printf("q, t, tubeid: %f %f %d \n",wcsimrootcherenkovdigihit->GetQ(),
                     wcsimrootcherenkovdigihit->GetT(),wcsimrootcherenkovdigihit->GetTubeId());
          }
      } // End of loop over Cherenkov digihits
#ifdef VERBOSE
      cout<<"calling BonsaiFit"<<endl;
#endif
      int vertexfound = bonsai->BonsaiFit( bsvertex, bsresult, bsgood, bsnsel, bsnhit, bsCAB, bsT, bsQ);
  //    BonsaiFit(float *vert, float *result, float *maxlike, int *nsel, int *nhit, int *cab, float *t, float *q)
  //        * int * nhit         : INPUT: Number of hits
  //        * int * cab          : INPUT: List of hit tube IDs. Length number of hits
  //        * float * t          : INPUT: List of hit times. Length number of hits
  //        * float * q          : INPUT: List of hit charges. Length number of hits
  //        * float * vert[4]    : OUTPUT: Reconstructed vertex (x,y,z,t)
  //        * float * result[5]  : OUTPUT: result[6] : Reconstructed direction = 
  //                                (zenith, azimuth, axial, cos(coneangle), ellipticity, goodness)
  //        * float * maxlike    : OUTPUT: [3] (goodness of time fit, goodness of position fit, maxq(?))
  //        * int * nsel         : OUTPUT: Number of selected hits
#ifdef VERBOSE
      (vertexfound) ? cout<<"fit successful"<<endl : cout<<"fit failed"<<endl;
#endif
      (vertexfound) ? numsuccessfulfits++ : numfailedfits++;
      if(true||vertexfound){  // record failed events so we can try to see why they failed
        //cout<<"vertex found"<<endl;
        // True vertex - don't fill until here in case bonsai failed
        (vertexfound) ? hTrueVtx_X->Fill(True_Vertex_X) : 0;
        (vertexfound) ? hTrueVtx_Y->Fill(True_Vertex_Y) : 0;
        (vertexfound) ? hTrueVtx_Z->Fill(True_Vertex_Z) : 0;
        (vertexfound) ? hTrueVtx_T->Fill(True_Vertex_T) : 0;
        (vertexfound) ? hTrueVtx_R->Fill(True_Vertex_R) : 0;
        
        // Reconstructed vertex
        Reco_Vertex_X = (vertexfound) ? bsvertex[0] : 0;
        Reco_Vertex_Y = (vertexfound) ? bsvertex[1] : 0;
        Reco_Vertex_Z = (vertexfound) ? bsvertex[2] : 0;
        Reco_Vertex_Vec=TVector3(Reco_Vertex_X,Reco_Vertex_Y,Reco_Vertex_Z);
        Reco_Vertex_T = (vertexfound) ? bsvertex[3]-950. : 0;
        Reco_Vertex_R=sqrt(pow(Reco_Vertex_X, 2) + pow(Reco_Vertex_Z, 2));
        (vertexfound) ? hRecoVtx_X->Fill(Reco_Vertex_X) : 0;
        (vertexfound) ? hRecoVtx_Y->Fill(Reco_Vertex_Y) : 0;
        (vertexfound) ? hRecoVtx_Z->Fill(Reco_Vertex_Z) : 0;
        (vertexfound) ? hRecoVtx_T->Fill(Reco_Vertex_T) : 0;
        (vertexfound) ? hRecoVtx_R->Fill(Reco_Vertex_R) : 0;
        vtxgoodness = (vertexfound) ? bsgood[2] : 0;
        (vertexfound) ? hVertexGoodness->Fill(vtxgoodness) : 0;
        
        // Vertex error
        Vertex_Diff_X = (vertexfound) ? True_Vertex_X-Reco_Vertex_X : 0;
        Vertex_Diff_Y = (vertexfound) ? True_Vertex_Y-Reco_Vertex_Y : 0;
        Vertex_Diff_Z = (vertexfound) ? True_Vertex_Z-Reco_Vertex_Z : 0;
        Vertex_Diff_T = (vertexfound) ? True_Vertex_T-Reco_Vertex_T : 0;
        Vertex_Diff_Vec=TVector3(Vertex_Diff_X,Vertex_Diff_Y,Vertex_Diff_Z);
        Vertex_Diff_Perp = (vertexfound) ? Vertex_Diff_Vec.Perp(True_Dir_Vec) : 0;
        Vertex_Diff_Par = (vertexfound) ? Vertex_Diff_Vec.Dot(True_Dir_Vec.Unit()) : 0;
        Vertex_Diff_Mag = (vertexfound) ? TMath::Sqrt(pow(Vertex_Diff_X,2)+pow(Vertex_Diff_Y,2)+pow(Vertex_Diff_Z,2)) : 0;
        (vertexfound) ? hVtxDiff_X->Fill(Vertex_Diff_X) : 0;
        (vertexfound) ? hVtxDiff_Y->Fill(Vertex_Diff_Y) : 0;
        (vertexfound) ? hVtxDiff_Z->Fill(Vertex_Diff_Z) : 0;
        (vertexfound) ? hVtxDiff_T->Fill(Vertex_Diff_T) : 0;
        (vertexfound) ? hVtxDiff_Par->Fill(Vertex_Diff_Par) : 0;
        (vertexfound) ? hVtxDiff_Perp->Fill(Vertex_Diff_Perp) : 0;
        (vertexfound) ? hVtxDiff_Tot->Fill(Vertex_Diff_Mag) : 0;

        // True direction
        (vertexfound) ? hTrueDir_X->Fill(True_Dir_X) : 0;
        (vertexfound) ? hTrueDir_Y->Fill(True_Dir_Y) : 0;
        (vertexfound) ? hTrueDir_Z->Fill(True_Dir_Z) : 0;
        
        // Reconstructed direction
        Reco_Dir_Zenith = (vertexfound) ? bsresult[0] : 0;
        Reco_Dir_Azimuth = (vertexfound) ? bsresult[1] : 0;
        Reco_Dir_Axial = (vertexfound) ? bsresult[2] : 0;
        Reco_Dir_CosOpAng = (vertexfound) ? bsresult[3] : 0;
        Reco_Dir_Ellip = (vertexfound) ? bsresult[4] : 0;
        dirgoodness = (vertexfound) ? bsresult[5] : 0;
        Reco_Dir_Vec = TVector3(0.,0.,1.);
        Reco_Dir_Vec.RotateY(Reco_Dir_Zenith);
        Reco_Dir_Vec.RotateZ(Reco_Dir_Azimuth);
        Reco_Dir_X = (vertexfound) ? Reco_Dir_Vec.X() : 0;
        Reco_Dir_Y = (vertexfound) ? Reco_Dir_Vec.Y() : 0;
        Reco_Dir_Z = (vertexfound) ? Reco_Dir_Vec.Z() : 0;
        (vertexfound) ? hRecoDir_X->Fill(Reco_Dir_X) : 0;
        (vertexfound) ? hRecoDir_Y->Fill(Reco_Dir_Y) : 0;
        (vertexfound) ? hRecoDir_Z->Fill(Reco_Dir_Z) : 0;
        (vertexfound) ? hDirectionGoodness->Fill(dirgoodness) : 0;
        
        // Direction error
        Dir_Diff_X = (vertexfound) ? True_Dir_X-Reco_Dir_X : 0;
        Dir_Diff_Y = (vertexfound) ? True_Dir_Y-Reco_Dir_Y : 0;
        Dir_Diff_Z = (vertexfound) ? True_Dir_Z-Reco_Dir_Z : 0;
        Dir_Diff_AngMag = (vertexfound) ? TMath::ACos((Reco_Dir_Vec.Unit()).Dot(True_Dir_Vec.Unit())) : 0;
        (vertexfound) ? hDirDiff_X->Fill(Dir_Diff_X) : 0;
        (vertexfound) ? hDirDiff_Y->Fill(Dir_Diff_Y) : 0;
        (vertexfound) ? hDirDiff_Z->Fill(Dir_Diff_Z) : 0;
        (vertexfound) ? hDirDiff_AngMag->Fill(Dir_Diff_AngMag) : 0;
        
        numdigitsthisevent=ncherenkovdigihits;
        totaltankcharge=wcsimrootevent->GetSumQ();
        (vertexfound) ? hNumDigits->Fill(ncherenkovdigihits) : 0;
        (vertexfound) ? hTotalCharge->Fill(totaltankcharge) : 0;
        
        treeout->Fill();
      } else { /*cout<<"vertex not found"<<endl;*/}
      
      // old code to calculate the SNR
      /*
      //cout<<"event "<<ev<<" had "<<ncherenkovdigihits<<" digits"<<endl;
      int numnoiseonlydigits=0, numpartialnoisedigits=0;
      int numdigitizedpesinthistrigger=0;
      // Loop over digits again
      for (int i=0;i<ncherenkovdigihits;i++){
        WCSimRootCherenkovDigiHit *thedigihit = 
          (WCSimRootCherenkovDigiHit*)wcsimrootevent->GetCherenkovDigiHits()->At(i);
        int thedigitstube = thedigihit->GetTubeId();
        Double_t thedigitscharge = thedigihit->GetQ();
        std::vector<int> photonids=thedigihit->GetPhotonIds();
        numdigitizedpesinthistrigger+=photonids.size();
        hNumHitsPerDigit->Fill(photonids.size());
        bool allnoise=true, somenoise=false;
        for(int counter=0; counter<photonids.size(); counter++){
          int thephotonsid=photonids.at(counter);
          WCSimRootCherenkovHitTime *thehittimeobject = 
            (WCSimRootCherenkovHitTime*)timeArray->At(thephotonsid);
          Int_t thephotonsparenttrackid=666;
          if(thehittimeobject) thephotonsparenttrackid = thehittimeobject->GetParentID();
          if(thephotonsparenttrackid==-1){ somenoise=true; } else { allnoise=false; }
          if(thephotonsparenttrackid==666){cout<<"digit "<<i<<", photon "<<counter<<" had no hittime"<<endl;}
        }
        if(somenoise&&!allnoise){
          //cout<<"DIGIT "<<i<<" had some noise hits *****************************"<<endl;
          numpartialnoisedigits++;
        } else if(allnoise){
          //cout<<"***************************************** DIGIT "<<i<<" ALL NOISE"<<endl;
          hNoiseCharge->Fill(thedigitscharge);
          numnoiseonlydigits++;
        } else if (!somenoise){
          //cout<<"no noise found in digit "<<i<<endl;
          hPhotonCharge->Fill(thedigitscharge);
        }
        hTotalCharge->Fill(thedigitscharge);
      } // End of loop over digi hits
      cout<<"number of pes in all digits in this trigger was "<<numdigitizedpesinthistrigger<<endl;
      //cout<<"of "<<ncherenkovdigihits<<" digits, "<<numpartialnoisedigits<<" digits had some noise hits, and "<<numnoiseonlydigits<<" digits were all noise"<<endl;
      //hNumNoiseOnlyDigits->Fill(numnoiseonlydigits);
      //hNumPartialNoiseDigits->Fill(numpartialnoisedigits);
      */
      
    } // End of loop over trigger (subevent)
    
    // reinitialize super event between loops.
#ifdef VERBOSE
    cout<<"reinitializing for next subtrigger"<<endl;
#endif
    wcsimrootsuperevent->ReInitialize();
    
  } // End of loop over events
  //  TCanvas c1("c1"); 
  
  cout<<"managed "<<numsuccessfulfits<<" successful fits but failed "<<numfailedfits<<" fits ("
      <<(numsuccessfulfits/(numsuccessfulfits+numfailedfits))*100<<"% success rate) and skipped"
      <<numskippedevents<<" events as no primary particle found"<<endl;
  
  Int_t bini=1;
  while(hVtxDiff_Tot->Integral(0,bini)<(hVtxDiff_Tot->GetEntries()*0.68)) bini++;
  Double_t one_sigma_res=hVtxDiff_Tot->GetBinLowEdge(bini);
  cout<<"one_sigma_res = "<<one_sigma_res<<endl;
  
  for(double runningsum=0, binno=0; binno<hVtxDiff_Tot->GetNbinsX(); binno++){
    runningsum+=hVtxDiff_Tot->GetBinContent(binno);
    hVtxDiff_Intg->SetBinContent(binno,runningsum);
  }
  double maxpercent=100.;
  hVtxDiff_Intg->Scale(maxpercent/hVtxDiff_Intg->Integral());
  
  for(double runningsum=0, binno=0; binno<hDirDiff_AngMag->GetNbinsX(); binno++){
    runningsum+=hDirDiff_AngMag->GetBinContent(binno);
    hDirDiff_Intg->SetBinContent(binno,runningsum);
  }
  hDirDiff_Intg->Scale(maxpercent/hDirDiff_Intg->Integral());
  
  // Set histogram titles, colours and draw.
  float win_scale = 0.75;
  int n_wide(2);
  int n_high(2);
  
  hNumDigits->GetXaxis()->SetTitle("Num Digits in Event");
  hNumDigits->GetYaxis()->SetTitle("Num Events");
  
  // Positions - true and reconstructed distributions
  hTrueVtx_X->SetLineColor(kRed);
  hTrueVtx_X->GetXaxis()->SetTitle("Vertex Position X [cm]");
  hTrueVtx_X->GetYaxis()->SetTitle("Num Events");
  hTrueVtx_Y->SetLineColor(kRed);
  hTrueVtx_Y->GetXaxis()->SetTitle("Vertex Position Y [cm]");
  hTrueVtx_Y->GetYaxis()->SetTitle("Num Events");
  hTrueVtx_Z->SetLineColor(kRed);
  hTrueVtx_Z->GetXaxis()->SetTitle("Vertex Position Z [cm]");
  hTrueVtx_Z->GetYaxis()->SetTitle("Num Events");
  hTrueVtx_R->SetLineColor(kRed);
  hTrueVtx_R->GetXaxis()->SetTitle("Vertex Radius [cm]");
  hTrueVtx_R->GetYaxis()->SetTitle("Num Events");
  hTrueVtx_T->SetLineColor(kRed);
  hTrueVtx_T->GetXaxis()->SetTitle("Vertex Position T [ns]");
  hTrueVtx_T->GetYaxis()->SetTitle("Num Events");
  hRecoVtx_X->SetLineColor(kBlue);
  hRecoVtx_X->GetXaxis()->SetTitle("Vertex Position X [cm]");
  hRecoVtx_X->GetYaxis()->SetTitle("Num Events");
  hRecoVtx_Y->SetLineColor(kBlue);
  hRecoVtx_Y->GetXaxis()->SetTitle("Vertex Position Y [cm]");
  hRecoVtx_Y->GetYaxis()->SetTitle("Num Events");
  hRecoVtx_Z->SetLineColor(kBlue);
  hRecoVtx_Z->GetXaxis()->SetTitle("Vertex Position Z [cm]");
  hRecoVtx_Z->GetYaxis()->SetTitle("Num Events");
  hRecoVtx_T->SetLineColor(kBlue);
  hRecoVtx_T->GetXaxis()->SetTitle("Vertex Time [ns]");
  hRecoVtx_T->GetYaxis()->SetTitle("Num Events");
  hRecoVtx_R->SetLineColor(kBlue);
  hRecoVtx_R->GetXaxis()->SetTitle("Vertex Radius [cm]");
  hRecoVtx_R->GetYaxis()->SetTitle("Num Events");
  // Draw the plots
  TLegend *leg1=nullptr; TCanvas *c1=nullptr, *c2=nullptr, *c3=nullptr, *c4=nullptr;
  if(DrawPlots){
    TCanvas* c1 = new TCanvas("c1", "First canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
    c1->Draw();
    c1->Divide(2,2);
    c1->cd(1); hTrueVtx_X->Draw();
    c1->cd(2); hTrueVtx_Y->Draw();
    c1->cd(3); hTrueVtx_Z->Draw();
    c1->cd(4); hTrueVtx_R->Draw();
    c1->cd(1); hRecoVtx_X->Draw("same");
    c1->cd(2); hRecoVtx_Y->Draw("same");
    c1->cd(3); hRecoVtx_Z->Draw("same");
    c1->cd(4); hRecoVtx_R->Draw("same");
    leg1 = new TLegend(0.2033426,0.8192982,0.45961,0.9210526,NULL,"brNDC");
    leg1->SetBorderSize(0);
    leg1->SetTextFont(62);
    leg1->SetTextSize(0.04093567);
    leg1->AddEntry(hTrueVtx_X,"True","l");
    leg1->AddEntry(hRecoVtx_X,"Reconstructed","l");
    c1->cd(1); leg1->Draw();
  }
  
  // Reconstructed position errors
  hVtxDiff_X->GetXaxis()->SetTitle("Reconstructed X Position Error [cm]");
  hVtxDiff_X->GetYaxis()->SetTitle("Num Events");
  hVtxDiff_Y->GetXaxis()->SetTitle("Reconstructed Y Position Error [cm]");
  hVtxDiff_Y->GetYaxis()->SetTitle("Num Events");
  hVtxDiff_Z->GetXaxis()->SetTitle("Reconstructed Z Position Error [cm]");
  hVtxDiff_Z->GetYaxis()->SetTitle("Num Events");
  hVtxDiff_T->GetXaxis()->SetTitle("Reconstructed T Error [ns]");
  hVtxDiff_T->GetYaxis()->SetTitle("Num Events");
  hVertexGoodness->GetXaxis()->SetTitle("Position Goodness of Fit");
  hVertexGoodness->GetYaxis()->SetTitle("Num Events");
  // Draw the plots
  if(DrawPlots){
    TCanvas* c2 = new TCanvas("c2", "Second canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
    c2->Divide(2,2);
    c2->cd(1); hVtxDiff_X->Draw();
    c2->cd(2); hVtxDiff_Y->Draw();
    c2->cd(3); hVtxDiff_Z->Draw();
    c2->cd(4); hVertexGoodness->Draw();
  }
  
  // Directions - Reconstructed and true distributions
  hTrueDir_X->SetLineColor(kRed);
  hTrueDir_X->GetXaxis()->SetTitle("Particle X Direction");
  hTrueDir_X->GetYaxis()->SetTitle("Num Events");
  hTrueDir_Y->SetLineColor(kRed);
  hTrueDir_Y->GetXaxis()->SetTitle("Particle Y Direction");
  hTrueDir_Y->GetYaxis()->SetTitle("Num Events");
  hTrueDir_Z->SetLineColor(kRed);
  hTrueDir_Z->GetXaxis()->SetTitle("Particle Z Direction");
  hTrueDir_Z->GetYaxis()->SetTitle("Num Events");
  hRecoDir_X->SetLineColor(kBlue);
  hRecoDir_X->GetXaxis()->SetTitle("Particle X Direction");
  hRecoDir_X->GetYaxis()->SetTitle("Num Events");
  hRecoDir_Y->SetLineColor(kBlue);
  hRecoDir_Y->GetXaxis()->SetTitle("Particle Y Direction");
  hRecoDir_Y->GetYaxis()->SetTitle("Num Events");
  hRecoDir_Z->SetLineColor(kBlue);
  hRecoDir_Z->GetXaxis()->SetTitle("Particle Z Direction");
  hRecoDir_Z->GetYaxis()->SetTitle("Num Events");
  // Draw the plots
  if(DrawPlots){
    TCanvas* c3 = new TCanvas("c3", "Third canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
    c3->Divide(2,2); 
    c3->cd(1); hTrueDir_X->Draw();
    c3->cd(2); hTrueDir_Y->Draw();
    c3->cd(3); hTrueDir_Z->Draw();
    c3->cd(1); hRecoDir_X->Draw("same");
    c3->cd(2); hRecoDir_Y->Draw("same");
    c3->cd(3); hRecoDir_Z->Draw("same");
    c3->cd(1); leg1->Draw();
  }
  
  // Reconstructed direction errors
  hDirDiff_X->GetXaxis()->SetTitle("Reconstructed X Direction Error");
  hDirDiff_X->GetYaxis()->SetTitle("Num Events");
  hDirDiff_Y->GetXaxis()->SetTitle("Reconstructed Y Direction Error");
  hDirDiff_Y->GetYaxis()->SetTitle("Num Events");
  hDirDiff_Z->GetXaxis()->SetTitle("Reconstructed Z Direction Error");
  hDirDiff_Z->GetYaxis()->SetTitle("Num Events");
  hDirDiff_AngMag->GetXaxis()->SetTitle("Direction Error Angular Diff");
  hDirDiff_AngMag->GetYaxis()->SetTitle("Num Events");
  hDirectionGoodness->GetXaxis()->SetTitle("Direction Goodness of Fit");
  hDirectionGoodness->GetYaxis()->SetTitle("Num Events");
  // Draw the plots
  if(DrawPlots){
    TCanvas* c4 = new TCanvas("c4", "Fourth canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
    c4->Divide(2,2); 
    c4->cd(1); hDirDiff_X->Draw();
    c4->cd(2); hDirDiff_Y->Draw();
    c4->cd(3); hDirDiff_Z->Draw();
    c4->cd(4); hDirDiff_AngMag->Draw();
  }
  
  // Write histograms to file
  fileout->cd();
  hTrueVtx_X->Write();
  hTrueVtx_Y->Write();
  hTrueVtx_Z->Write();
  hTrueVtx_T->Write();
  hTrueVtx_R->Write();
  hRecoVtx_X->Write();
  hRecoVtx_Y->Write();
  hRecoVtx_Z->Write();
  hRecoVtx_T->Write();
  hRecoVtx_R->Write();
  hVtxDiff_X->Write();
  hVtxDiff_Y->Write();
  hVtxDiff_Z->Write();
  hVtxDiff_T->Write();
  hVtxDiff_Tot->Write();
  hVtxDiff_Par->Write();
  hVtxDiff_Perp->Write();
  hVtxDiff_Intg->Write();
  hVertexGoodness->Write();
  hTrueDir_X->Write();
  hTrueDir_Y->Write();
  hTrueDir_Z->Write();
  hRecoDir_X->Write();
  hRecoDir_Y->Write();
  hRecoDir_Z->Write();
  hDirDiff_X->Write();
  hDirDiff_Y->Write();
  hDirDiff_Z->Write();
  hDirDiff_AngMag->Write();
  hDirDiff_Intg->Write();
  hDirectionGoodness->Write();
  
  hNumDigits->Write();
//  hNumHitsPerDigit->Write();
//  hNumNoiseOnlyDigits->Write();
//  hNumPartialNoiseDigits->Write();
//  hNoiseCharge->Write();
//  hPhotonCharge->Write();
  hTotalCharge->Write();
  treeout->Write();
  
  fileout->Close();
  delete fileout;
  if(DrawPlots){
    if(c1) delete c1; c1=0;
    if(c2) delete c2; c2=0;
    if(c3) delete c3; c3=0;
    if(c4) delete c4; c4=0;
    if(leg1) delete leg1; leg1=0;
  }
  
  return 0;
}
