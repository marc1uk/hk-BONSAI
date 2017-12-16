/* vim:set noexpandtab tabstop=2 wrap */
//C++
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <sys/types.h> // for stat() test to see if file or folder
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
//ROOT
#include <TH1F.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TChain.h>
#include <TStopwatch.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TH1.h>
#include <TRint.h>
#include <TColor.h>
#include <TString.h>
//#include "../gitver/Bonsai_v0/bonsai/pmt_geometry.h"
//#include "../gitver/Bonsai_v0/bonsai/likelihood.h"
//#include "../gitver/Bonsai_v0/bonsai/goodness.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
//WCSim
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
//BONSAI
#include "WCSimBonsai.hh"
// WCSim Reconstruction
#include "../../../WCSim/gitver/root_work/MRDspecs.hh"
#endif

#ifndef VERBOSE
//#define VERBOSE 1
#endif


template <typename T>
std::vector<T> valsinrange(typename std::vector<T>::iterator rangebegin, typename std::vector<T>::iterator rangeend, T lowerlimit, T upperlimit){
  auto lowerit = std::lower_bound(rangebegin, rangeend, lowerlimit);
  auto upperit = std::upper_bound(rangebegin, rangeend, upperlimit);
  return typename std::vector<T> (lowerit, upperit);
}

// anonymous namespace to declare constants used within this translation unit.
namespace {
  constexpr float LIGHT_SPEED = 21.58333; // speed of light in water, value from https://github.com/hyperk/hk-BONSAI/blob/d9b227dad26fb63f2bfe80f60f7f58b5a703250a/bonsai/hits.h#L5
  constexpr int MAXDIGITS=500; // maximum number of digits we allow, due to use of preallocated vector sizes.
  constexpr int nWorkingPMTs = -1; // number of working PMTs If negative, = nPMTs from WCSimRootGeom
  constexpr float darkRate = 0;//4.2/1000000; // dark noise rate (per ns) of the PMT (dummy value, based on 8.4kHz for B&L PMT)
  constexpr float lambdaEff = 100*100; // scattering length in cm (dummy value, based on Design Report II.2.E.1)
  constexpr float INTER_PMT_SEPARATION_CM = 106.f; // distance to neighboring PMTs is 106cm azimuthally, slightly less for ring separation, except at the end rings, which are further
  constexpr float effCoverages[] = {0.4, 0.4, 0.4, 0.4, 0.4068, 0.4244, 0.4968, 0.5956, 0.67};
  // from MC: effective photodetector coverage at an incidence angle of theta = 5, 15, ..., 85 degree
  constexpr double triggeroffset = 950.0; // WCSimRootTrigger::offset. TODO put in WCSimRootOptions
  constexpr bool ZERO_DIGIT_TS = false;   // whether to subtract off triggertime + triggeroffset before
                                          // handling digits so that they're relative to reference time 0.
  
  constexpr int track_pdg = 11; // Which track PDG to extract for Truth comparisons.
  // the highest energy particle of this kind will be selected.
  constexpr bool primaries_only = true; // the matched track must be a primary.
  constexpr bool isANNIE = true; // enable some modifications specific for ANNIE (e.g. z is the BEAM AXIS)
  constexpr bool DrawPlots=false; // draw plots at the end of the event loop.
  
  // function to calculate the reconstructed energy from the number of 'effective' hits:
  // this is approximately linear, except at very low energies, but the fit form and values
  // FIXME FIXME should be calculated your detector! FIXME FIXME
  // note a constexpr has many restrictions: it is possible this qualifier may need to be lifted if you change the form.
  constexpr float EnergyFromnEff(float nEff){
    return ((25.00 + 0.138*(nEff-189.8))*0.378);
//    float a[5]= {0.82, 0.13, -1.11*pow(10, -4), 1.25*pow(10, -6), -3.42*pow(10, -9)};
//    if (nEff<189.8) {
//      for (int n=0;n<5;n++) {
//        eRec += a[n]*pow(nEff, n);
//      }
//    } else {
//      eRec=(25.00 + 0.138*(nEff-189.8))*0.378;
//    }
  }
  
  constexpr int BOGUS_INT = std::numeric_limits<int>::max();
  constexpr float BOGUS_FLOAT = std::numeric_limits<float>::max();
  constexpr double BOGUS_DOUBLE = std::numeric_limits<double>::max();
}

void GetTankExitPoint(const TVector3 RecoVtx_Vec, const TVector3 RecoDir_Vec, TVector3 &projectedinterceptpoint, double &tracklengthintank, WCSimRootGeom* geo);

// XXX NOTE: XXX ALL RESULTS ARE WITH THE TANK CENTRED AT (0,0,0) EXCEPT FOR: projectedmrdentrypoint and tankexitpoint

//const double miploss = 2.5;
//const double ecorroffset = 204.105;
//const double ecorrgrad = -2.52276;
//const TF1 energylossvslength("lin","pol1(0)");
////TF1 energylossvslength("lin","pol1(0)");
////energylossvslength.SetParameters(204.105,-2.52276)

// Simple example of reading a generated Root file
int WCSim_Bonsai_Ana(const char *filedirectory=NULL, bool verbose=false)
{

//////////////////////////////////////////////////////////////////////////////////////////////
// Load Libraries and Header
//////////////////////////////////////////////////////////////////////////////////////////////
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
  
//////////////////////////////////////////////////////////////////////////////////////////////
// Create Output File and Tree
//////////////////////////////////////////////////////////////////////////////////////////////
  //TODO: make the tree entry into a class, so we can easily friend the trees and grab the entry
  // as an object. Maybe use MakeClass to implement this.
  TFile* fileout = new TFile("bonsaiout.root","RECREATE");
  TTree* treeout = new TTree("bonsaitree","Bonsai Reconstruction Tree");
  // Info to locate the event
  std::string wcsimfilepath;
  Int_t run_id;
  Int_t event_id;
  Int_t subtrigger;
  TBranch* bWCSimFile = treeout->Branch("WCSim_File",&wcsimfilepath);
  TBranch* bRunId = treeout->Branch("RunId",&run_id);
  TBranch* bEventId = treeout->Branch("EventId",&event_id);
  TBranch* bSubtriggerId = treeout->Branch("SubtriggerId",&subtrigger);
  // info about the event
  double numdigitsthisevent; // ncherenkovdigihits
  TBranch* bNumDigits = treeout->Branch("NumDigits",&numdigitsthisevent);
  double totaltankcharge;
  TBranch* bSumQ = treeout->Branch("TotalTankCharge",&totaltankcharge);
  // Vertex
  double True_Vertex_X, True_Vertex_Y, True_Vertex_Z, True_Vertex_T, True_Vertex_R;
  TBranch* bTrueVtx_X = treeout->Branch("TrueVtx_X",&True_Vertex_X);
  TBranch* bTrueVtx_Y = treeout->Branch("TrueVtx_Y",&True_Vertex_Y);
  TBranch* bTrueVtx_Z = treeout->Branch("TrueVtx_Z",&True_Vertex_Z);
  TBranch* bTrueVtx_R = treeout->Branch("TrueVtx_T",&True_Vertex_T);
  TBranch* bTrueVtx_T = treeout->Branch("TrueVtx_R",&True_Vertex_R);
  TVector3 True_Vertex_Vec;
  TBranch* bTrueVtx_Vec = treeout->Branch("TrueVtx_Vec",&True_Vertex_Vec);
  int vertexfound;
  TBranch* bVertex_Found = treeout->Branch("Vertex_Found",&vertexfound);
  bool vertexintank;
  TBranch* bVertexInTank = treeout->Branch("Vertex_In_Tank",&vertexintank);
  double Reco_Vertex_X, Reco_Vertex_Y, Reco_Vertex_Z, Reco_Vertex_T, Reco_Vertex_R;
  TBranch* bRecoVtx_X = treeout->Branch("RecoVtx_X",&Reco_Vertex_X);
  TBranch* bRecoVtx_Y = treeout->Branch("RecoVtx_Y",&Reco_Vertex_Y);
  TBranch* bRecoVtx_Z = treeout->Branch("RecoVtx_Z",&Reco_Vertex_Z);
  TBranch* bRecoVtx_T = treeout->Branch("RecoVtx_T",&Reco_Vertex_T);
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
  double vtxgoodness;
  TBranch* bVtx_Goodness = treeout->Branch("Vtx_Goodness",&vtxgoodness);
  // Direction
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
  double RecoDir_Error;
  TBranch* bRecoDir_Error = treeout->Branch("RecoDir_Error",&RecoDir_Error);
  double Dir_Diff_X, Dir_Diff_Y, Dir_Diff_Z, Dir_Diff_AngMag;
  TBranch* bDirDiff_X = treeout->Branch("DirDiff_X",&Dir_Diff_X);
  TBranch* bDirDiff_Y = treeout->Branch("DirDiff_Y",&Dir_Diff_Y);
  TBranch* bDirDiff_Z = treeout->Branch("DirDiff_Z",&Dir_Diff_Z);
  TBranch* bDirDiff_AngMag = treeout->Branch("DirDiff_AngMag",&Dir_Diff_AngMag);
  double dirgoodness;
  TBranch* bDir_Goodness = treeout->Branch("Dir_Goodness",&dirgoodness);
  double fittingtime;
  TBranch* bFittingTime = treeout->Branch("Fitting_Time",&fittingtime);
  std::vector<double> numusedhits;
  TBranch* bNumUsedHits = treeout->Branch("NHitsUsedForFit",&numusedhits);
  TVector3 tankexitpoint;
  TBranch* bTankExitPoint = treeout->Branch("Tank_Exit_Point",&tankexitpoint);
  double tracklengthintank;
  TBranch* bTrackLengthInTank = treeout->Branch("TrackLengthInTank",&tracklengthintank);
  double tankenergyloss;
  TBranch* bTankEnergyLoss = treeout->Branch("TankEnergyLoss",&tankenergyloss);
  double tankenergylosserror;
  TBranch* bTankEnergyLossError = treeout->Branch("TankEnergyLossError",&tankenergylosserror);
  bool mrdintercept;
  TBranch* bMrdIntercept = treeout->Branch("Intercepts_MRD",&mrdintercept);
  TVector3 projectedmrdentrypoint;
  TBranch* bMrdEntryPoint = treeout->Branch("Projected_MRDEntry",&projectedmrdentrypoint);
  Int_t MrdTrackId;
  TBranch* bMrdTrackId = treeout->Branch("MrdTrackId",&MrdTrackId);
  // Energy
  double True_Particle_E_ini, True_Particle_Mom_ini, True_Particle_E_fin=BOGUS_DOUBLE, True_Particle_Mom_fin=BOGUS_DOUBLE, Reco_Particle_E_loss;
  TBranch* bTrueE = treeout->Branch("TrueEini",&True_Particle_E_ini);
  TBranch* bTrueMom = treeout->Branch("TrueMomMagIni",&True_Particle_Mom_ini);
  TBranch* bTrueEfin = treeout->Branch("TrueEfin",&True_Particle_E_fin);
  TBranch* bTrueMomfin = treeout->Branch("TrueMomMagFin",&True_Particle_Mom_fin);
  TBranch* bRecoE = treeout->Branch("RecoEloss",&Reco_Particle_E_loss); // based on digits / effective hits
  double N_Effective_Hits;
  TBranch* bNeffectiveHits = treeout->Branch("NeffectiveHits",&N_Effective_Hits);
  
//////////////////////////////////////////////////////////////////////////////////////////////
// Create Output Histograms
//////////////////////////////////////////////////////////////////////////////////////////////
  
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
  
  TH1F *recEnergy = new TH1F("Reconstructed Energy Loss", "Reconstructed Energy Loss;Energy Loss (MeV); NEvents", 50, 0, 100);
  TH1F *recEnergyDiff = new TH1F("Reconstructed Energy Loss Error", "Reconstructed Energy Loss  Error;Energy Loss Error Ratio [MeV]); NEvents", 50, 0, 100);
  TH1F *recEnergyDiffRatio = new TH1F("Reconstructed Energy Loss Error Ratio", "Reconstructed Energy Loss  Error;Energy Loss Error Ratio [no units]); NEvents", 50, 0, 100);
  
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
  
//////////////////////////////////////////////////////////////////////////////////////////////
// Load File and Tree
//////////////////////////////////////////////////////////////////////////////////////////////
  gROOT->cd();
  // Declare input paths, files, trees...
  TFile* wcsimfile=0;
  TString wcsimfilename, wcsimfilestring;
  TTree* wcsimT=0;
  
  // first check if the path we've been given is a file or folder
  bool isdir, addsubfolders=true;
  struct stat s;
  if(stat(filedirectory,&s)==0){
    if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory?? how does this work?
      isdir=true;  //it's a directory
    } else if(s.st_mode & S_IFREG){ // mask to check if it's a file??
      isdir=false; //it's a file
    } else {
      assert(false&&"Check input path: stat says it's neither file nor directory..?");
    }
  } else {
    //assert(false&&"stat failed on input path! Is it valid?"); // error
    // errors could be because this is a file pattern: e.g. wcsim_0.4*.root Treat as a file.
    isdir=false;
  }
  
  // TChain for wcsim files
  TChain* c =  new TChain("wcsimT");
  TString chainpattern;
  if(isdir){
    chainpattern = TString::Format("%s/wcsim_0.1000.root",filedirectory);
    //TString chainpattern = TString::Format("%s/wcsim_10MeV_iso_e_wDN_Multidigit_16-04-17.root",filedirectory);
  } else {
    chainpattern = filedirectory;
  }
  
  cout<<"loading TChain entries from "<<chainpattern<<endl;
  c->Add(chainpattern);
  Int_t nevent = c->GetEntries();
  cout<<"loaded "<<nevent<<" entries in the chain"<<endl;
  if(nevent<1){ return -9; }
  
  // wcsimT
  WCSimRootEvent* wcsimrootsuperevent=0, *m=0, *v=0;
  TBranch* bp=0, *mp=0, *vp=0;
  WCSimRootTrigger* wcsimrootevent=0, *atrigm=0, *atrigv=0;
  WCSimRootGeom* geo = 0; 
  int numpmts, nWorkingPMTsToUse;
  
  cout<<"loading first wcsimT tree from "<<chainpattern<<" tchain"<<endl;
  c->LoadTree(0);
  Int_t treeNumber = -1;
  wcsimT = c->GetTree();
  Int_t thistreesentries = wcsimT->GetEntries();
  cout<<thistreesentries<<" entries in the first tree"<<endl;
  
  WCSimBonsai* bonsai = new WCSimBonsai();
  
//////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN LOOP OVER EVENTS
//////////////////////////////////////////////////////////////////////////////////////////////
  // Now loop over events
  //nevent=10;
  int numsuccessfulfits=0, numfailedfits=0, numskippedevents=0;
#ifdef VERBOSE
  cout<<nevent<<" events to analyze"<<endl;
#endif
  for (int ev=0; ev<nevent; ev++){
  
#ifdef VERBOSE
    cout<<"loading entry "<<ev<<endl;
#endif
    Long64_t localEntry = c->LoadTree(ev);
    if( localEntry<0){ cout<<"end of tchain"<<endl; break; }
    Int_t nextTreeNumber = c->GetTreeNumber();
    if(treeNumber!=nextTreeNumber){
      cout<<"new tree: "<<nextTreeNumber<<endl;
      // this means we've switched file - need to load the new branch addresses (and potentially other trees)
      // first pull out the new file name
      wcsimT = c->GetTree();
      wcsimfile = wcsimT->GetCurrentFile();
      wcsimfilename=wcsimfile->GetName();
      wcsimfilestring=wcsimfilename+std::string(filedirectory);
      Int_t thistreesentries = wcsimT->GetEntries();
      cout<<"wcsimT has "<<thistreesentries<<" entries in this file"<<endl;
  
      // load the geometry tree and grab the geometry if we haven't already
      if(geo==0){
//        TString geofileloc = TString::Format("%s/../wcsim_wdirt_17-06-17/wcsim_0.1000.root",filedirectory);
//        TString geofileloc = "/home/marc/LinuxSystemFiles/Bonsai/gitver/Bonsai_v0/wcsim_0.1001.root";
//        TFile* geofile = TFile::Open(geofileloc.Data());
        TFile* geofile = wcsimfile;
        if(geofile==nullptr) assert(false&&"geofile doesn't exist!");
        TTree* geotree = (TTree*)geofile->Get("wcsimGeoT");
        if(geotree==0){ cout<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
        geotree->SetBranchAddress("wcsimrootgeom", &geo);
        if (geotree->GetEntries() == 0) { cout<<"geotree has no entries!"<<endl; exit(9); }
        geotree->GetEntry(0);
        bonsai->Init(geo, isANNIE);
        numpmts = geo->GetWCNumPMT();
        nWorkingPMTsToUse = (nWorkingPMTs<0) ? numpmts : nWorkingPMTs;
        // pmtidtocopynum = makecopynummap(geo); // turns out this is a 1:1 mapping after all!
        //MakePMTmap(geo, topcappositionmap, bottomcappositionmap, wallpositionmap);
      }
  
      // wcsim trigger classes
      wcsimT->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent, &bp);
      //wcsimT->SetBranchAddress("wcsimrootevent_mrd",&m, &mp);
      //wcsimT->SetBranchAddress("wcsimrootevent_facc",&v, &vp);
      bp->SetAutoDelete(kTRUE);
      //mp->SetAutoDelete(kTRUE);
      //vp->SetAutoDelete(kTRUE);
      if(bp==0/*||mp==0||vp==0*/){ cout<<"branches are zombies!"<<endl; }
  
      treeNumber=nextTreeNumber;
    }
  if(nextTreeNumber>0) verbose=true;
  if(treeNumber==0&&localEntry>3) continue;
   //////////////////////////////
   
    // Read the event from the tree into the WCSimRootEvent instance
    wcsimT->GetEntry(localEntry);
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
    
//////////////////////////////////////////////////////////////////////////////////////////////
// Get True Track Info
//////////////////////////////////////////////////////////////////////////////////////////////
    
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
      if( (primaries_only && nextrack->GetParenttype()==0 && nextrack->GetFlag()==0) &&
          (nextrack->GetIpnu()==track_pdg&&nextrack->GetE()>Emax) ){
        Emax=nextrack->GetE();
        mutrackindex=track;
      }
    }
#ifdef VERBOSE
    cout<<"mu found"<<endl;
#endif
    if(mutrackindex<0) {numskippedevents++; if(localEntry>10) assert(false&&"no mu"); continue;}// this event had no muons!
    cout<<"✩ ";
    if(verbose){
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0), 
         wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    } else {
      //std::cout<<"Event "<<ev<<std::endl;
    }
#ifdef VERBOSE
    cout<<"numtruehits ="<< wcsimrootevent->GetCherenkovHits()->GetEntries()<<endl;
    cout<<"numdigits = "<< wcsimrootevent->GetCherenkovDigiHits()->GetEntries()<<endl;
#endif
    
    nextrack = (WCSimRootTrack*)wcsimrootevent->GetTracks()->At(mutrackindex);
    True_Particle_E_ini = nextrack->GetE();
    True_Particle_Mom_ini = nextrack->GetP();
#if FILE_VERSION>3
    True_Particle_E_fin = nextrack->GetEndE();
    True_Particle_Mom_fin = nextrack->GetEndP();
#endif
    if(isANNIE){
      True_Vertex_X=nextrack->GetStart(0)-geo->GetWCOffset(0);
      True_Vertex_Y=nextrack->GetStart(1)-geo->GetWCOffset(1);
      True_Vertex_Z=nextrack->GetStart(2)-geo->GetWCOffset(2);
      //cout<<"chosen nextrack starts at: ("<<nextrack->GetStart(0)<<","
      //    <<nextrack->GetStart(1)<<", "<<nextrack->GetStart(2)<<")"<<endl
      //    <<"adding offset: ("<<True_Vertex_X<<", "<<True_Vertex_X<<", "<<True_Vertex_X<<")"<<endl;
      True_Vertex_R=sqrt(pow(True_Vertex_X, 2) + pow(True_Vertex_Z, 2));
      True_Vertex_T=nextrack->GetTime();
      True_Vertex_Vec=TVector3(True_Vertex_X,True_Vertex_Y,True_Vertex_Z);
      
      True_Dir_X=nextrack->GetDir(0);
      True_Dir_Y=nextrack->GetDir(1);
      True_Dir_Z=nextrack->GetDir(2);
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

//////////////////////////////////////////////////////////////////////////////////////////////
// Perform Bonsai Fit
//////////////////////////////////////////////////////////////////////////////////////////////
    
    // Loop over sub events and perform vertex fit
    // ===========================================
    // static arrays for the fit
    float bsVertex[4],bsResult[6];
    float bsGood[3];
    int *bsnhit; //nsel (SLE)
    int bsNsel[2]; //nsel (SLE)
    std::vector<float> bsT(MAXDIGITS,0), bsQ(MAXDIGITS,0), distance(MAXDIGITS,0), PMTX(MAXDIGITS,0), PMTY(MAXDIGITS,0), PMTZ(MAXDIGITS,0), tCorrected(MAXDIGITS,0);
    std::vector<int> n50Array(MAXDIGITS,0), bsCAB(MAXDIGITS,0);

    if(verbose) std::cout << "DIGITIZED HITS:" << std::endl;
    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
    for (int index = 0 ; index < 1 /*wcsimrootsuperevent->GetNumberOfEvents()*/; index++){
    // TODO: for tree entries to tie up with WCSim etc, we need one entry per WCSim event
    // disable analysis of delayed triggers until the class supports multiple triggers
    // in one tree entry
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
      if(verbose) std::cout << "Sub event number = " << index << "\n";
      WCSimRootEventHeader* header = wcsimrootevent->GetHeader();
      int triggertime = header->GetDate();
#ifdef VERBOSE
      if(ZERO_DIGIT_TS) cout<<"adding offset to digits: + "<<triggertime<<" - "<<triggeroffset<<endl;
#endif
      
      int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
//      cout<<"numtruehits = "<<wcsimrootevent->GetCherenkovHits()->GetEntries()
//          <<", numdigits = "<<wcsimrootevent->GetCherenkovDigiHits()->GetEntries()<<endl;
      if(ncherenkovdigihits==0) continue;
      if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      bsnhit = & ncherenkovdigihits;
#ifdef VERBOSE
      cout<<"filling array of Ts and Qs"<<endl;
#endif
      for (int i=0;i<ncherenkovdigihits;i++){
          // Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
          TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
          
          WCSimRootCherenkovDigiHit *cherenkovdigihit = 
            dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
          
          if(ZERO_DIGIT_TS){
            bsT.at(i) = cherenkovdigihit->GetT() + triggertime - triggeroffset;
          } else {
            bsT.at(i) = cherenkovdigihit->GetT();
          }
          bsQ.at(i) = cherenkovdigihit->GetQ();
          bsCAB.at(i) = cherenkovdigihit->GetTubeId(); // N.B. SHOULD start from 1. used in hits::hits.
          
          WCSimRootPMT pmt = geo->GetPMT(cherenkovdigihit->GetTubeId()-1);
          if(not isANNIE){
            PMTX.at(i) = pmt.GetPosition(0)-geo->GetWCOffset(0);
            PMTY.at(i) = pmt.GetPosition(1)-geo->GetWCOffset(1);
            PMTZ.at(i) = pmt.GetPosition(2)-geo->GetWCOffset(2);
          } else {
            PMTX.at(i) = pmt.GetPosition(0)-geo->GetWCOffset(0);
            PMTY.at(i) = pmt.GetPosition(2)-geo->GetWCOffset(2);
            PMTZ.at(i) = pmt.GetPosition(1)-geo->GetWCOffset(1);
          }
          
          if(verbose){
            if ( i < 10 ) // Only print first XX=10 tubes
              printf("q, t, tubeid: %f %f %d \n",cherenkovdigihit->GetQ(),
                     cherenkovdigihit->GetT(),cherenkovdigihit->GetTubeId());
          }
      } // End of loop over Cherenkov digihits
#ifdef VERBOSE
      cout<<"calling BonsaiFit"<<endl;
#endif
      TStopwatch* ms = new TStopwatch();
      ms->Start();
//      cout<<"calling BonsaiFit with "<<bsCAB.size()<<" bsCAB: {";
//      for(auto atube : bsCAB) cout<<atube<<", ";
//      cout<<"}"<<endl<<"bsT: {";
//      for(auto atime : bsT) cout<<atime<<", ";
//      cout<<"}"<<endl<<"bsQ: {";
//      for(auto acharge : bsQ) cout<<acharge<<", ";
//      cout<<"}"<<endl;
      vertexfound = bonsai->BonsaiFit( bsVertex, bsResult, bsGood, bsNsel, bsnhit, bsCAB.data(), bsT.data(), bsQ.data());
//      cout<<"bonsaifit result is:"<<endl<<"bsVertex: {";
//      for(int i=0; i<4; i++) cout<<bsVertex[i]<<", ";
//      cout<<"}"<<endl<<"bsResult: {";
//      for(int i=0; i<5; i++) cout<<bsResult[i]<<", ";
//      cout<<"}"<<endl<<"bsGood: {";
//      for(int i=0; i<3; i++) cout<<bsGood[i]<<", ";
//      cout<<"}"<<endl<<"bsNsel: {";
//      for(int i=0; i<2; i++) cout<<bsNsel[i]<<", ";
//      cout<<"}"<<endl<<"bsnhit: "<<*bsnhit<<endl;
      ms->Stop();
      fittingtime = ms->CpuTime();
      numusedhits.clear();
      for(int i=0; i<2; i++) numusedhits.push_back(bsNsel[i]); // what does this really store?*****
#ifdef VERBOSE
      cout << " Vertex Fitting Time :  Real = " << ms->RealTime() << " ; CPU = " << ms->CpuTime() << "\n";
#endif
  //    BonsaiFit(float *vert, float *result, float *maxlike, int *nsel, int *nhit, int *cab, float *t, float *q)
  //        * int * nhit         : INPUT: Number of hits
  //        * int * cab          : INPUT: List of hit tube IDs. Length number of hits
  //        * float * t          : INPUT: List of hit times. Length number of hits
  //        * float * q          : INPUT: List of hit charges. Length number of hits
  //        * float * vert[4]    : OUTPUT: Reconstructed vertex (x,y,z,t)
  //        * float * result[5]  : OUTPUT: result[6] : Reconstructed direction = 
  //                                (zenith, azimuth, axial, cos(coneangle), ellipticity, goodness)
  //        * float * maxlike    : OUTPUT: [3] (goodness of time fit, goodness of position fit, maxq(?))
  //        * int * nsel         : OUTPUT: Number of selected hits XXX BUT IT'S AN ARRAY OF 2 VALS?****
#ifdef VERBOSE
      (vertexfound) ? cout<<"fit successful"<<endl : cout<<"fit failed"<<endl;
#endif
      (vertexfound) ? numsuccessfulfits++ : numfailedfits++;
      
//////////////////////////////////////////////////////////////////////////////////////////////
// Retrieve Bonsai Vertex Information
//////////////////////////////////////////////////////////////////////////////////////////////
      
      if(true||vertexfound){  // record failed events so we can try to see why they failed
        //cout<<"vertex found"<<endl;
        
        // Reconstructed vertex
        Reco_Vertex_X = (vertexfound) ? bsVertex[0] : 0;
        if(isANNIE){
          Reco_Vertex_Y = (vertexfound) ? bsVertex[2] : 0;
          Reco_Vertex_Z = (vertexfound) ? bsVertex[1] : 0;
        } else {
          Reco_Vertex_Y = (vertexfound) ? bsVertex[1] : 0;
          Reco_Vertex_Z = (vertexfound) ? bsVertex[2] : 0;
        }
        Reco_Vertex_Vec=TVector3(Reco_Vertex_X,Reco_Vertex_Y,Reco_Vertex_Z);
        if(not ZERO_DIGIT_TS){
          Reco_Vertex_T = (vertexfound) ? bsVertex[3]-triggeroffset : 0;
           // FIXME should this also account for triggertime? ^^^^
        } else {
          Reco_Vertex_T = (vertexfound) ? bsVertex[3] : 0;
        }
        Reco_Vertex_R=sqrt(pow(Reco_Vertex_X, 2) + pow(Reco_Vertex_Z, 2));
        vtxgoodness = (vertexfound) ? bsGood[2] : 0;
        
        // Vertex error - only set if we have a valid reconstructed vertex
        Vertex_Diff_X = (vertexfound) ? True_Vertex_X-Reco_Vertex_X : 0;
        Vertex_Diff_Y = (vertexfound) ? True_Vertex_Y-Reco_Vertex_Y : 0;
        Vertex_Diff_Z = (vertexfound) ? True_Vertex_Z-Reco_Vertex_Z : 0;
        Vertex_Diff_T = (vertexfound) ? True_Vertex_T-Reco_Vertex_T : 0;
        Vertex_Diff_Vec=TVector3(Vertex_Diff_X,Vertex_Diff_Y,Vertex_Diff_Z);
        Vertex_Diff_Perp = (vertexfound) ? TMath::Abs(Vertex_Diff_Vec.Perp(True_Dir_Vec)) : 0;
        Vertex_Diff_Par = (vertexfound) ? TMath::Abs(Vertex_Diff_Vec.Dot(True_Dir_Vec.Unit())) : 0;
        Vertex_Diff_Mag = (vertexfound) ? TMath::Sqrt(pow(Vertex_Diff_X,2)+pow(Vertex_Diff_Y,2)+pow(Vertex_Diff_Z,2)) : 0;
        
        // Reconstructed direction
        Reco_Dir_Zenith = (vertexfound) ? bsResult[0] : 0;
        Reco_Dir_Azimuth = (vertexfound) ? bsResult[1] : 0;
        Reco_Dir_Axial = (vertexfound) ? bsResult[2] : 0;
        Reco_Dir_CosOpAng = (vertexfound) ? bsResult[3] : 0;
        Reco_Dir_Ellip = (vertexfound) ? bsResult[4] : 0;
        dirgoodness = (vertexfound) ? bsResult[5] : 0;
        // Conversion from Euler angles to cartesian unit vector
        Reco_Dir_Vec = TVector3(0.,1.,0.);
        Reco_Dir_Vec.SetMagThetaPhi(1.,Reco_Dir_Zenith,Reco_Dir_Azimuth);
        Reco_Dir_X = (vertexfound) ? Reco_Dir_Vec.X() : 0;
        Reco_Dir_Y = (vertexfound) ? Reco_Dir_Vec.Z() : 0;
        Reco_Dir_Z = (vertexfound) ? Reco_Dir_Vec.Y() : 0;
        RecoDir_Error = (vertexfound) ? TMath::Pi()/4. : 0.; // TODO figure out how to calculate this !!
        
        // Direction error
        Dir_Diff_X = (vertexfound) ? True_Dir_X-Reco_Dir_X : 0;
        Dir_Diff_Y = (vertexfound) ? True_Dir_Y-Reco_Dir_Y : 0;
        Dir_Diff_Z = (vertexfound) ? True_Dir_Z-Reco_Dir_Z : 0;
        Dir_Diff_AngMag = (vertexfound) ? TMath::ACos((Reco_Dir_Vec.Unit()).Dot(True_Dir_Vec.Unit())) : 0;
        
        // Other event info
        vertexintank = (vertexfound) ? 
          ((Reco_Vertex_R<tank_radius)&&(abs(Reco_Vertex_Y)<tank_halfheight)) : false;
        numdigitsthisevent=ncherenkovdigihits;
        totaltankcharge=wcsimrootevent->GetSumQ();
        (vertexfound) ? hNumDigits->Fill(ncherenkovdigihits) : 0;
        (vertexfound) ? hTotalCharge->Fill(totaltankcharge) : 0;
        
        
//////////////////////////////////////////////////////////////////////////////////////////////
// Energy estimation for this trigger
//////////////////////////////////////////////////////////////////////////////////////////////
      // For a detailed description of the energy estimation formulas in SK-IV, see ch. 4.3 of
      // http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf

      for (int i=0;i<ncherenkovdigihits;i++) {// Loop through all WCSimRootCherenkovDigiHits in this trigger
        // get distance of hit (=PMT position) to reconstructed vertex (bsVertex[i])
        distance.at(i) = sqrt(pow((PMTX.at(i)-bsVertex[0]), 2) + pow((PMTY.at(i)-bsVertex[1]), 2) + pow((PMTZ.at(i)-bsVertex[2]), 2));

        // substract time-of-flight from measured arrival time bsT[i] --> tCorrected.at(i)
        tCorrected.at(i) = bsT.at(i) - (distance.at(i)/LIGHT_SPEED);
      }
      auto tCorrectedEnd = tCorrected.begin()+ncherenkovdigihits;

      // sort tCorrected array into ascending order
      std::sort(tCorrected.begin(),tCorrectedEnd);

      // look for the 50 ns interval with the maximum total number of hits --> start time: tMin
      float tMin;
      int n50Max=0;
      int iValue=0;

      for (int i=0;i<ncherenkovdigihits;i++) { // loop through tCorrected array: take each element as tMin and find 50ns window
        tMin = tCorrected.at(i);
        auto hitsintimerange = valsinrange(tCorrected.begin(),tCorrectedEnd,tMin,tMin+50.f);
        if(hitsintimerange.size() > n50Max){
          n50Max = hitsintimerange.size();
          iValue = i;
        }
      }
      tMin = tCorrected.at(iValue);

      // find the number of hits in the 100 ns interval corresponding to the maximal 50 ns window
      auto hitsintimerange = valsinrange(tCorrected.begin(),tCorrectedEnd,tMin,tMin+100.f);
      int n100Max = hitsintimerange.size();
//      cout<<"digits in the range "<<tMin<<" - "<<(tMin+100.f)<<" are: (";
//      for(auto aval : hitsintimerange) cout<<aval<<", ";
//      cout<<")"<<endl;
//      cout<<"all digits are: (";
//      for(int digitcount=0; digitcount<ncherenkovdigihits; digitcount++)
//        cout<<tCorrected.at(digitcount)<<", ";
//      cout<<")"<<endl;

      int j=0;
      std::vector<float> distance50(n50Max);
      std::vector<int> tubeID(n50Max);
      // fill arrays of distance from vertex in cm and tubeID for each hit in maximal interval
      for (int i=0;i<ncherenkovdigihits;i++) { // loop each hit in ncherenkovdigihits
        if (tMin <=tCorrected.at(i) && tCorrected.at(i) < tMin + 50.f) {
          distance50.at(j) = sqrt(pow((PMTX.at(i)-bsVertex[0]), 2) + pow((PMTY.at(i)-bsVertex[1]), 2) + pow((PMTZ.at(i)-bsVertex[2]), 2));
          tubeID.at(j) = bsCAB.at(i);
          j++;
        }
      } // end of loop over hits

      if (verbose) {
        std::cout << "tMin: " << tMin << "\n";
        std::cout << "n50Max: " << n50Max << "\n";
        std::cout << "n100Max: " << n100Max << "\n";
      }

      float nEff = 0; // effective number of hits
      float occupancy;
      for (int i=0;i<n50Max;i++) { // loop over hits in 50 ns interval and calculate nEff
        WCSimRootPMT pmt = geo->GetPMT(tubeID.at(i)-1);
//        cout<<"get pmt at "<<i<<endl;
        float x,y,z;
        if(not isANNIE){
          x = pmt.GetPosition(0)-geo->GetWCOffset(0);
          y = pmt.GetPosition(1)-geo->GetWCOffset(1);
          z = pmt.GetPosition(2)-geo->GetWCOffset(2);
        } else {
          x = pmt.GetPosition(0)-geo->GetWCOffset(0);
          y = pmt.GetPosition(2)-geo->GetWCOffset(2);
          z = pmt.GetPosition(1)-geo->GetWCOffset(1);
        }
//        cout<<"pos is ("<<x<<", "<<y<<", "<<z<<")"<<endl;

        // calculate occupancy to correct for multiple hits on a single PMT
        // In a 3x3 grid around PMT 'tubeID', what proportion of PMTs has seen a hit?
        // TODO: Treat PMTs at the edge (that have fewer neighbors) differently!
        int nearbyHits = 0;
        WCSimRootPMT otherPMT;
        for (int j=0; j<n50Max; j++) { // loop through all hit PMTs, count number of hits in nearby PMTs
//          cout<<"get other pmt at "<<j<<"="<<(tubeID.at(j)-1)<<endl;
          otherPMT = geo->GetPMT(tubeID.at(j)-1);
          double distancetoother;
          if(not isANNIE){
            distancetoother = sqrt(pow(x - (otherPMT.GetPosition(0)-geo->GetWCOffset(0)), 2) + pow(y - (otherPMT.GetPosition(1)-geo->GetWCOffset(1)), 2) + pow(z - (otherPMT.GetPosition(2)-geo->GetWCOffset(2)), 2));
          } else {
            distancetoother = sqrt(pow(x - (otherPMT.GetPosition(0)-geo->GetWCOffset(0)), 2) + pow(y - (otherPMT.GetPosition(2)-geo->GetWCOffset(2)), 2) + pow(z - (otherPMT.GetPosition(1)-geo->GetWCOffset(1)), 2));
          }
//          cout<<"pos is ("<<otherPMT.GetPosition(0)-geo->GetWCOffset(0)
//            <<", "<<otherPMT.GetPosition(2)-geo->GetWCOffset(2)
//            <<", "<<otherPMT.GetPosition(1)-geo->GetWCOffset(1)
//            <<"), distance is "<<distancetoother<<endl;
          if (distancetoother < INTER_PMT_SEPARATION_CM) {
            nearbyHits++;
          }
        }

        double ratio = float(nearbyHits) / 9.;
        if (ratio < 1.) {
          occupancy= log(1. / (1.-ratio)) / ratio; // from Poisson statistics
        } else {
          occupancy= 3.0;
        }

        // correct for delayed hits (e.g. due to scattering)
        
        float lateHits = (n100Max - n50Max - (nWorkingPMTsToUse * darkRate * 50.f)) / float(n50Max);

        // substract dark noise hits
        float darkNoise = (nWorkingPMTsToUse * darkRate * 50.f) / float(n50Max);

        // calculate effective coverage to correct for photoCathodeCoverage
        // this depends on angle of incidence, see Fig. 4.5 (left) of http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf
        // TODO: take into account azimuthal dependence as show in Fig. 4.5 (right); currently assumes phi=0
//        cout<<"occupancy="<<occupancy<<", lateHits="<<lateHits<<", darkNoise="<<darkNoise;
        float dotProduct;
        if(not isANNIE){
          dotProduct = pmt.GetOrientation(0)*(bsVertex[0] - x) + pmt.GetOrientation(1)*(bsVertex[1] - y) + pmt.GetOrientation(2)*(bsVertex[2] - z);
        } else {
          dotProduct = pmt.GetOrientation(0)*(bsVertex[0] - x) + pmt.GetOrientation(2)*(bsVertex[1] - y) + pmt.GetOrientation(1)*(bsVertex[2] - z);
        }
//        cout<<", dotProduct="<<dotProduct;
        float theta = acos( dotProduct / distance50.at(i)) * 180.f / TMath::Pi();
        if (theta > 89.99) theta = 0.f; // we have apparently mis-reconstructed the vertex, so let's set ...
        if (theta < 0.f) theta = 0.f; // ... the coverage to the most likely value of 0.4 (i.e. theta < 40 degrees)
//        float phi = 0; // dummy value
//        cout<<", theta="<<theta;
        float photoCathodeCoverage = 1 / effCoverages[int(theta/10)];
//        cout<<"photoCathodeCoverage="<<photoCathodeCoverage<<endl;

        // correct for scattering in water
        float waterTransparency = exp(distance50.at(i) / lambdaEff);

        float nEffHit = (occupancy + lateHits - darkNoise) * photoCathodeCoverage * waterTransparency;
        nEff += nEffHit;

        if (verbose) {
          std::cout << "\n*** event #" << ev << ", PMT hit #" << i << " ***************************************\n";
          std::cout << "occupancy (ratio of hits in 3x3 grid): " << occupancy << " (" << ratio << ")\n";
          std::cout << "lateHits:  " << lateHits << "\n";
          std::cout << "darkNoise: " << darkNoise << "\n";
          std::cout << "photoCathodeCoverage: " << photoCathodeCoverage << "\n";
          std::cout << "waterTransparency: " << waterTransparency << "\n";
          std::cout << "nEff for this 1 hit: " << nEffHit << "\n";
        }
      } // end of loop over hits in 50 ns interval


      nEff *= numpmts / float(nWorkingPMTsToUse); // correct for dead PMTs; convert nWorkingPMTs to float because integer division is inaccurate

      float eRec = EnergyFromnEff(nEff);
      recEnergy->Fill(eRec);
      recEnergyDiff->Fill(TMath::Abs(True_Particle_E_fin-True_Particle_E_ini-eRec));
      recEnergyDiffRatio->Fill(TMath::Abs(True_Particle_E_fin-True_Particle_E_ini-eRec)/(True_Particle_E_fin-True_Particle_E_ini));
      
      Reco_Particle_E_loss = (double)eRec;
      N_Effective_Hits = (double)nEff;

      if (verbose) {
        std::cout << "Neff = " << nEff << std::endl;
        std::cout << "Reconstructed energy = " << eRec << std::endl;
      }
        
//////////////////////////////////////////////////////////////////////////////////////////////
// Fill Histograms and Tree
//////////////////////////////////////////////////////////////////////////////////////////////
        
        // True vertex - always fill this info
        hTrueVtx_X->Fill(True_Vertex_X);
        hTrueVtx_Y->Fill(True_Vertex_Y);
        hTrueVtx_Z->Fill(True_Vertex_Z);
        hTrueVtx_T->Fill(True_Vertex_T);
        hTrueVtx_R->Fill(True_Vertex_R);
        // True direction
        hTrueDir_X->Fill(True_Dir_X);
        hTrueDir_Y->Fill(True_Dir_Y);
        hTrueDir_Z->Fill(True_Dir_Z);
        
        // Reconstructed Vertex - only fill if found
        (vertexfound) ? hRecoVtx_X->Fill(Reco_Vertex_X) : 0;
        (vertexfound) ? hRecoVtx_Y->Fill(Reco_Vertex_Y) : 0;
        (vertexfound) ? hRecoVtx_Z->Fill(Reco_Vertex_Z) : 0;
        (vertexfound) ? hRecoVtx_T->Fill(Reco_Vertex_T) : 0;
        (vertexfound) ? hRecoVtx_R->Fill(Reco_Vertex_R) : 0;
        (vertexfound) ? hVertexGoodness->Fill(vtxgoodness) : 0;
        
        // Vertex error
        (vertexfound) ? hVtxDiff_X->Fill(Vertex_Diff_X) : 0;
        (vertexfound) ? hVtxDiff_Y->Fill(Vertex_Diff_Y) : 0;
        (vertexfound) ? hVtxDiff_Z->Fill(Vertex_Diff_Z) : 0;
        (vertexfound) ? hVtxDiff_T->Fill(Vertex_Diff_T) : 0;
        (vertexfound) ? hVtxDiff_Par->Fill(Vertex_Diff_Par) : 0;
        (vertexfound) ? hVtxDiff_Perp->Fill(Vertex_Diff_Perp) : 0;
        (vertexfound) ? hVtxDiff_Tot->Fill(Vertex_Diff_Mag) : 0;
        
        // Reconstructed Direction
        (vertexfound) ? hRecoDir_X->Fill(Reco_Dir_X) : 0;
        (vertexfound) ? hRecoDir_Y->Fill(Reco_Dir_Y) : 0;
        (vertexfound) ? hRecoDir_Z->Fill(Reco_Dir_Z) : 0;
        (vertexfound) ? hDirectionGoodness->Fill(dirgoodness) : 0;
        
        // Direction Error
        (vertexfound) ? hDirDiff_X->Fill(Dir_Diff_X) : 0;
        (vertexfound) ? hDirDiff_Y->Fill(Dir_Diff_Y) : 0;
        (vertexfound) ? hDirDiff_Z->Fill(Dir_Diff_Z) : 0;
        (vertexfound) ? hDirDiff_AngMag->Fill(Dir_Diff_AngMag) : 0;
        
        // Project to tank exit RecoVtx_X
        if(vertexfound) GetTankExitPoint(Reco_Vertex_Vec, Reco_Dir_Vec, tankexitpoint, tracklengthintank, geo);
        //TODO tankexitxbounds, tankexitybounds;
        //     how to calculate these, as we have no error on fit vtx/dir?
        
        // TODO: above calculates track length in tank. Derive energy loss from this, or from charge deposition
        tankenergyloss=tracklengthintank*2.5; // nominally 2.5MeV/cm
        tankenergylosserror=50.;  //TODO: base on length error, rate error
                                  // OR redo this based on charge deposition fit fit error
        if(vertexfound){
          double mrdentryx = Reco_Vertex_X+((Reco_Dir_X/Reco_Dir_Z)*(MRD_start-(Reco_Vertex_Z+geo->GetWCOffset(2))));
          double mrdentryy = Reco_Vertex_Y+((Reco_Dir_Y/Reco_Dir_Z)*(MRD_start-(Reco_Vertex_Z+geo->GetWCOffset(2))));
          projectedmrdentrypoint=TVector3(mrdentryx,mrdentryy,MRD_start);
          mrdintercept = ((abs(mrdentryx)<MRD_width)&&(abs(mrdentryy)<MRD_height));
        } else {
          projectedmrdentrypoint=TVector3(0.,0.,0.);
          mrdintercept=false;
        }
        
        // not useful here but will be when this tree entry becomes a class object
        MrdTrackId = -1;
        
        // Set any remaining variables
        wcsimfilepath = wcsimfilestring;
        run_id = wcsimrootevent->GetHeader()->GetRun();
        event_id = ev;
        Int_t subtrigger= index;
        
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
    
//////////////////////////////////////////////////////////////////////////////////////////////
// Calculate Aggregate Info
//////////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<"managed "<<numsuccessfulfits<<" successful fits but failed "<<numfailedfits<<" fits ("
      <<(((double)numsuccessfulfits/((double)numsuccessfulfits+(double)numfailedfits))*100.)
      <<"% success rate) and skipped "<<numskippedevents<<" events with no true primary particle found"<<endl;
  
  Int_t bini=1;
  while(hVtxDiff_Tot->Integral(0,bini)<(hVtxDiff_Tot->GetEntries()*0.68)) bini++;
  Double_t one_sigma_res=hVtxDiff_Tot->GetBinLowEdge(bini);
  cout<<"one_sigma_res = "<<one_sigma_res<<endl;
  
  double totalentries = hVtxDiff_Tot->GetEntries();
  for(double runningsum=0, binno=1; binno<hVtxDiff_Tot->GetNbinsX()+1; binno++){
    runningsum+=hVtxDiff_Tot->GetBinContent(binno);
    hVtxDiff_Intg->SetBinContent(binno,(runningsum/totalentries));
  }
  
  for(double runningsum=0, binno=1; binno<hDirDiff_AngMag->GetNbinsX()+1; binno++){
    runningsum+=hDirDiff_AngMag->GetBinContent(binno);
    hDirDiff_Intg->SetBinContent(binno,(runningsum/totalentries));
  }
  
//////////////////////////////////////////////////////////////////////////////////////////////
// Set Histogram Properties
//////////////////////////////////////////////////////////////////////////////////////////////
  
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
  
//////////////////////////////////////////////////////////////////////////////////////////////
// Write Histograms to File
//////////////////////////////////////////////////////////////////////////////////////////////
  
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
  recEnergy->Write();
  recEnergyDiff->Write();
  recEnergyDiffRatio->Write();
  
  hNumDigits->Write();
//  hNumHitsPerDigit->Write();
//  hNumNoiseOnlyDigits->Write();
//  hNumPartialNoiseDigits->Write();
//  hNoiseCharge->Write();
//  hPhotonCharge->Write();
  hTotalCharge->Write();
  
//////////////////////////////////////////////////////////////////////////////////////////////
// Write Tree to File and Close
//////////////////////////////////////////////////////////////////////////////////////////////
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

void GetTankExitPoint(const TVector3 RecoVtx_Vec, const TVector3 RecoDir_Vec, TVector3 &projectedinterceptpoint, double &tracklengthintank, WCSimRootGeom* geo){
// NOTE: TANK IS CENTRED ON (0,0,0) FOR RECONSTRUCTED VERTEX
  // ----------------------------------------------------------------------------------------------
  // calculate the track length in water
  // ----------------------------------------------------------------------------------------------
  // to calculate track length _in water_ find distance from start vertex to the point
  // where it intercepts the tank. if this length > total track length, use total track length
  // otherwise use this length. 

  if( (abs(RecoVtx_Vec.X())>tank_radius) || (abs(RecoVtx_Vec.Z())>tank_radius) ||
      (sqrt(pow(RecoVtx_Vec.X(),2.)+pow(RecoVtx_Vec.Z(),2.))>tank_radius) ||
      (abs(RecoVtx_Vec.Y())>tank_halfheight) ){
      projectedinterceptpoint = TVector3(BOGUS_FLOAT,BOGUS_FLOAT,BOGUS_FLOAT);
      tracklengthintank = -1.;
      return; // not gonna intercept the tank if it's already out it.
  }
  double trackgradx = RecoDir_Vec.X()/RecoDir_Vec.Z();
  double trackgrady = RecoDir_Vec.Y()/RecoDir_Vec.Z();
  double xatziszero = RecoVtx_Vec.X() - ((RecoVtx_Vec.Z())*(trackgradx));
#ifdef VERBOSE
  cout<<"RecoVec is ("<<RecoVtx_Vec.X()<<", "<<RecoVtx_Vec.Y()<<", "<<RecoVtx_Vec.Z()<<")"<<endl;
  cout<<"RecoDir is ("<<RecoDir_Vec.X()<<", "<<RecoDir_Vec.Y()<<", "<<RecoDir_Vec.Z()
      <<"), gradients are "<<trackgradx<<", "<<trackgrady<<endl;
  cout<<"projected to centre of tank, X value is "<<xatziszero<<endl;
#endif
  double firstterm = -trackgradx*xatziszero;
  double thirdterm = 1+TMath::Power(trackgradx,2.);
  double secondterm = (TMath::Power(tank_radius,2.)*thirdterm) - TMath::Power(xatziszero,2.);
  if(secondterm<0){ // for safety
    cerr<<"failure of vertex reconstruction! Reconstructed track does not project to intercept tank radius!"<<endl;
    cerr<<"RecoVec is ("<<RecoVtx_Vec.X()<<", "<<RecoVtx_Vec.Y()<<", "<<RecoVtx_Vec.Z()<<")"<<endl;
    cerr<<"RecoDir is ("<<RecoDir_Vec.X()<<", "<<RecoDir_Vec.Y()<<", "<<RecoDir_Vec.Z()<<")"<<endl;
    cerr<<"Tank radius: "<<tank_radius<<", tank halfheight: "<<tank_halfheight<<endl;
    projectedinterceptpoint=TVector3(0.,0.,0.);
    tracklengthintank=0.;
    assert(false);
  }
  double solution1 = (firstterm + TMath::Sqrt(secondterm))/thirdterm;
  double solution2 = (firstterm - TMath::Sqrt(secondterm))/thirdterm;
  double tankendpointz;
  if(RecoDir_Vec.Z()>0){
    tankendpointz = solution1;  //forward going track
  } else {
    tankendpointz = solution2;  // backward going track
  }
  double tankendpointx = RecoVtx_Vec.X() + (tankendpointz-RecoVtx_Vec.Z())*trackgradx;
  // now check if the particle would have exited through one of the caps before reaching this radius
  double tankendpointy = RecoVtx_Vec.Y() + (tankendpointz-RecoVtx_Vec.Z())*trackgrady;
  double tankendpointx1=tankendpointx;
  double tankendpointz1=tankendpointz;
  double tankendpointy1=tankendpointy;
  double distanceoutsidecap =TMath::Abs(tankendpointy)-(tank_halfheight);
  if(distanceoutsidecap>0 && 
      TMath::Abs(RecoVtx_Vec.Y())<(tank_halfheight)){
    // this trajectory exits through the cap. Need to recalculate x, z exiting points...!
    if(RecoDir_Vec.Y()>0){
      tankendpointy = tank_halfheight;  // by definition of leaving through cap
    } else {
      tankendpointy = -tank_halfheight;
    }
    tankendpointz = RecoVtx_Vec.Z() + (tankendpointy-RecoVtx_Vec.Y())/trackgrady;
    tankendpointx = RecoVtx_Vec.X() + (tankendpointz-RecoVtx_Vec.Z())*trackgradx;
  } else {
    // this trajectory exited the tank by a side point; existing value is valid
  }
  // check if new reconstructed radius is outside the radial extent
  double distanceoutsidebarrel = sqrt(pow(tankendpointx,2)+pow(tankendpointz,2))-(tank_radius);
  if(distanceoutsidebarrel>0){
    // this also seems to have failed. We could just bail out:
//    cerr<<"failure of vertex reconstruction! Reconstructed track lies outside tank radius!"<<endl;
//    cerr<<"RecoVec is ("<<RecoVtx_Vec.X()<<", "<<RecoVtx_Vec.Y()<<", "<<RecoVtx_Vec.Z()<<")"<<endl;
//    cerr<<"RecoDir is ("<<RecoDir_Vec.X()<<", "<<RecoDir_Vec.Y()<<", "<<RecoDir_Vec.Z()<<")"<<endl;
//    projectedinterceptpoint=TVector3(0.,0.,0.);
//    tracklengthintank=0.;
//    return;
    // nah. just use whichever vertex is closest to the tank edges
    if(distanceoutsidebarrel>distanceoutsidecap){
    // go back to old points
    tankendpointx=tankendpointx1;
    tankendpointy=tankendpointy1;
    tankendpointz=tankendpointz1;
    }
  }
  
  // we're now able to determine muon track length in the tank:
  tracklengthintank = TMath::Sqrt(
    TMath::Power((tankendpointx-RecoVtx_Vec.X()),2)+
    TMath::Power((tankendpointy-RecoVtx_Vec.Y()),2)+
    TMath::Power((tankendpointz-RecoVtx_Vec.Z()),2) );

  Double_t maxtanktracklength = 
  TMath::Sqrt(TMath::Power(tank_radius*2.,2.)+TMath::Power(tank_halfheight*2.,2.));
  if(tracklengthintank > maxtanktracklength || TMath::IsNaN(tracklengthintank)){
    cerr<<"Reconstructed vertex is ("<<RecoVtx_Vec.X()<<", "<<RecoVtx_Vec.Y()<<", "<<RecoVtx_Vec.Z()
        <<"), tank endpoint is ("<<tankendpointx<<", "<<tankendpointy<<", "<<tankendpointz
        <<"), distanceoutsidecap="<<distanceoutsidecap
        <<", distanceoutsidebarrel="<<distanceoutsidebarrel<<endl;
    if(distanceoutsidecap>0){
      cerr<<"original radial exit vertex was ("<<tankendpointx1<<", "<<tankendpointy1<<", "<<tankendpointz1;
    }
    cerr<<"Track length "<<tracklengthintank<<", max possible length "<<maxtanktracklength<<"!"<<endl;
    cerr<<"Tank radius: "<<tank_radius<<", tank halfheight: "<<tank_halfheight<<endl;
    if(tracklengthintank > maxtanktracklength){
      cerr<<"track length is too long!"<<endl;
    } else {
      cerr<<"NaN RESULT FROM TRACK LENGTH IN TANK?!"<<endl;
    }
    assert(false);
  }

  // correct for tank y,z offset
  tankendpointz += geo->GetWCOffset(2);
  tankendpointy += geo->GetWCOffset(1);
  projectedinterceptpoint=TVector3(tankendpointx,tankendpointy,tankendpointz);
}
