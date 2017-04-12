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

// Simple example of reading a generated Root file
int sample_bonsai(const char *filename=NULL, bool verbose=false)
{
  // Clear global scope
  //gROOT->Reset();
  /*
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetTitleColor(1);
    gStyle->SetStatColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleSize(0.04);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPalette(1);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(.5);
    gStyle->SetTitleY(0.99);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetHatchesLineWidth(2);
    gStyle->SetLineWidth(1.5);
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetTitleSize(0.04,"X");
    gStyle->SetTitleSize(0.04,"Y");
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasBorderMode(0);
  */
  
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

  TFile* fileout = new TFile("bonsaiout.root","RECREATE");
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

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;

//  ===============================================================
//  ANNIE WCSim variables: tank -152<X<152, -212<Y<183, 15<Z<320 cm
//  ===============================================================
  const Float_t tank_start = 15.70;          // front face of the tank in cm
  const Float_t tank_radius = 152.4;         // tank radius in cm
  const Float_t tank_halfheight = 198.;      // tank half height in cm
  const Float_t tank_yoffset = -14.46;       // tank y offset in cm
  
  TH1F *hTrueVtx_X = new TH1F("Event True VTX_X", "Event True VTX_X", 200, -400, 400);
  TH1F *hTrueVtx_Y = new TH1F("Event True VTX_Y", "Event True VTX_Y", 200, -400, 400);
  TH1F *hTrueVtx_Z = new TH1F("Event True VTX_Z", "Event True VTX_Z", 200, -300, 300);
  //TH1F *hTrueVtx_T = new TH1F("Event True VTX_T", "Event True VTX_T", 200, -1500, 1500);
  TH1F *hRecoVtx_X = new TH1F("Reconstructed X", "Reconsructed X", 200, -400, 400);
  TH1F *hRecoVtx_Y = new TH1F("Reconstructed Y", "Reconsructed Y", 200, -400, 400);
  TH1F *hRecoVtx_Z = new TH1F("Reconstructed Z", "Reconsructed Z", 200, -300, 300);
  TH1F *hRecoVtx_T = new TH1F("Reconstructed T", "Reconsructed T", 200, -100, 100);
  TH1F *hvtxR = new TH1F("Reconstructed R", "Reconsructed R", 200, 0, 600);
  TH1F *bsgy = new TH1F("Bonsai Goodness", "Bonsai Goodness", 200, 0, 1);
  //TH1F *hRecoVtx_Y = new TH1F("Event VTX2", "Event VTX2", 200, -1500, 1500);
  //TH1F *hRecoVtx_Y = new TH1F("SumQ", "SumQ", 200, 0, 200);
  TH1F *hVtxDiff_X = new TH1F("Vertex Error X", "Vertex Error X", 200, -200, 200);
  TH1F *hVtxDiff_Y = new TH1F("Vertex Error Y", "Vertex Error Y", 200, -200, 200);
  TH1F *hVtxDiff_Z = new TH1F("Vertex Error Z", "Vertex Error Z", 200, -200, 200);
  
  TH1F *hNumDigits = new TH1F("Num PMT Digits", "Num PMT Digits", 200, 0, 200);

  // Now loop over events
  for (int ev=0; ev<nevent; ev++)
    {
      // Read the event from the tree into the WCSimRootEvent instance
      tree->GetEntry(ev);      
      wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
      if(verbose){
	printf("********************************************************");
	printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	       wcsimrootevent->GetHeader()->GetDate());
	printf("Mode %d\n", wcsimrootevent->GetMode());
	printf("Number of subevents %d\n",
	       wcsimrootsuperevent->GetNumberOfSubEvents());
	
	printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
	// use GetVtxs(1,X) not GetVtx(X) to take first secondary, not neutrino, in case we didn't have GENIE
	// if the event is in the tank, they will be the same
	printf("Vtx %f %f %f\n", wcsimrootevent->GetVtxs(1,0),
	       wcsimrootevent->GetVtxs(1,1),wcsimrootevent->GetVtxs(1,2));
      } else {
	std::cout<<"Event "<<ev<<std::endl;
      }
      Float_t True_Vertex_X=wcsimrootevent->GetVtx(0)*1000.;
      Float_t True_Vertex_Y=wcsimrootevent->GetVtx(1)*1000.;
      Float_t True_Vertex_Z=wcsimrootevent->GetVtx(2)*1000.;
/*
      if(std::to_string(True_Vertex_X).substr(0,5)=="-99.9"){ // we don't have neutrino vtx info, need to load tracks:
        int numtracks= wcsimrootevent->GetNtrack();
        // scan through the truth tracks, find the primary muon and pull vertex info from it
        for(int track=0; track<numtracks; track++){
          WCSimRootTrack* nextrack = (WCSimRootTrack*)wcsimrootevent->GetTracks()->At(track);
//           a WCSimRootTrack has methods: 
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
//          Float_t   GetTime()             trj->GetGlobalTime(); stopping(?) time of particle
//          Int_t     GetId()               wcsim trackid
//          
          if(nextrack->GetParenttype()==0){ 
            True_Vertex_X=nextrack->GetStart(0);
            True_Vertex_Y=nextrack->GetStart(1);
            True_Vertex_Z=nextrack->GetStart(2);
            break;
          }
        }
      }
*/
//    don't fill histograms with true vertex until we've checked if we managed a fit
      
      if(verbose){
	printf("Jmu %d\n", wcsimrootevent->GetJmu());
	printf("Npar %d\n", wcsimrootevent->GetNpar());
	printf("Ntrack %d\n", wcsimrootevent->GetNtrack());
      }
      // Now read the tracks in the event
      
      // Get the number of tracks
      int ntrack = wcsimrootevent->GetNtrack();
      if(verbose) printf("ntracks=%d\n",ntrack);
      
      int i;
      // Loop through elements in the TClonesArray of WCSimTracks
      for (i=0; i<ntrack; i++)
	{
	  TObject *element = (wcsimrootevent->GetTracks())->At(i);
	  
	  WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
	  
	  if(verbose){
	    printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	    printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
	    
	    for (int j=0; j<3; j++)
	      printf("Track dir: %d %f\n",j, wcsimroottrack->GetDir(j));
	  }
	  
	  
	}  // End of loop over tracks
      
      // Now look at the Cherenkov hits

      // Get the number of Cherenkov hits.
      // Note... this is *NOT* the number of photons that hit tubes.
      // It is the number of tubes hit with Cherenkov photons.
      // The number of digitized tubes will be smaller because of the threshold.
      // Each hit "raw" tube has several photon hits.  The times are recorded.
      // See http://nwg.phy.bnl.gov/DDRD/cgi-bin/private/ShowDocument?docid=245
      // for more information on the structure of the root file.
      //  
      // The following code prints out the hit times for the first 10 tubes and also
      // adds up the total pe.
      // 
      // For digitized info (one time/charge tube after a trigger) use
      // the digitized information.
      //
      
      int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
      int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits(); 
      
      hNumDigits->Fill(ncherenkovdigihits);
      if(verbose){
	printf("node id: %i\n", ev);
	printf("Ncherenkovhits %d\n",     ncherenkovhits);
	printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
	std::cout << "RAW HITS:" << std::endl;
      }

      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();
      
      int totalPe = 0;
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      for (i=0; i< ncherenkovhits; i++)
	{
	  TObject *Hit = (wcsimrootevent->GetCherenkovHits())->At(i);
	  WCSimRootCherenkovHit *wcsimrootcherenkovhit = 
	    dynamic_cast<WCSimRootCherenkovHit*>(Hit);
	  
	  int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	  int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	  int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	  WCSimRootPMT pmt   = geo->GetPMT(tubeNumber-1);
	  totalPe += peForTube;
	  
	  
	  if ( i < 10 ) // Only print first XX=10 tubes
	    {
	      if(verbose) printf("Total pe: %d times( ",peForTube);
	      for (int j = timeArrayIndex; j < timeArrayIndex + peForTube; j++)
		{
		  WCSimRootCherenkovHitTime * HitTime = 
		    dynamic_cast<WCSimRootCherenkovHitTime*>(timeArray->At(j));
		  
		  if(verbose) printf("%6.2f ", HitTime->GetTruetime() );	     
		}
	      if(verbose) std::cout << ")" << std::endl;
	    }
	  
	} // End of loop over Cherenkov hits
      if(verbose) std::cout << "Total Pe : " << totalPe << std::endl;
      
      // Look at digitized hit info
      
      // Get the number of digitized hits
      // Loop over sub events
      float bsT[500],bsQ[500];
      float bsvertex[4],bsresult[6];
      float bsgood[500];
      int bsCAB[500];
      int bsnhit[1]; //nsel (SLE)
      int bsnsel[2]; //nsel (SLE)
      if(verbose) std::cout << "DIGITIZED HITS:" << std::endl;
      for (int index = 0 ; index < /*wcsimrootsuperevent->GetNumberOfEvents()*/1; index++) 
	{
	  wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
	  if(verbose) std::cout << "Sub event number = " << index << "\n";
	  
	  int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
	  if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
	  
	  //for (i=0;i<(ncherenkovdigihits>4 ? 4 : ncherenkovdigihits);i++){
	  bsnhit[0] = ncherenkovdigihits;
	  for (i=0;i<ncherenkovdigihits;i++)
	    {
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
//  BonsaiFit(float *vert, float *result, float *maxlike, int *nsel, int *nhit, int *cab, float *t, float *q)
//  * int * nhit         : INPUT: Number of hits
//  * int * cab          : INPUT: List of hit tube IDs. Length number of hits
//  * float * t          : INPUT: List of hit times. Length number of hits
//  * float * q          : INPUT: List of hit charges. Length number of hits
//  * float * vert[4]    : OUTPUT: Reconstructed vertex (x,y,z,t)
//  * float * result[5]  : OUTPUT: result[5] : Reconstructed direction (x,y,z,t,ll0)
//                         t is ...? ll0 is the likelihood of that reconstructed direction.
//  * float * maxlike    : OUTPUT: [2] is goodness of position fit. [1] is goodness of time fit
//  * int * nsel         : OUTPUT: Number of selected hits
	  int vertexfound = bonsai->BonsaiFit( bsvertex, bsresult, bsgood, bsnsel, bsnhit, bsCAB, bsT, bsQ);
	  if(vertexfound){  // 
	    hRecoVtx_X->Fill(bsvertex[0]);
	    hRecoVtx_Y->Fill(bsvertex[1]);
	    hRecoVtx_Z->Fill(bsvertex[2]);
	    hRecoVtx_T->Fill(bsvertex[3]-950.);
	    hvtxR->Fill(sqrt(pow(bsvertex[0], 2) + pow(bsvertex[1], 2) + pow(bsvertex[2], 2)));
	    bsgy->Fill(bsgood[2]);
	    
	    hVtxDiff_X->Fill(True_Vertex_X-bsvertex[0]);
	    hVtxDiff_Y->Fill((True_Vertex_Y+tank_yoffset)-bsvertex[1]);
	    hVtxDiff_Z->Fill((True_Vertex_Z-tank_start-tank_radius)-bsvertex[2]);
	    
	    hTrueVtx_X->Fill(True_Vertex_X);
	    hTrueVtx_Y->Fill(True_Vertex_Y+tank_yoffset);
	    hTrueVtx_Z->Fill(True_Vertex_Z-tank_start-tank_radius);
	    //hTrueVtx_T->Fill(wcsimrootevent->???);  // not directly saved in WCSim
	    //hRecoVtx_Y->Fill(wcsimrootevent->GetSumQ());
	  }
	  
	} // End of loop over trigger
      
      // reinitialize super event between loops.
      wcsimrootsuperevent->ReInitialize();
      
    } // End of loop over events
  //  TCanvas c1("c1"); 
  float win_scale = 0.75;
  int n_wide(2);
  int n_high(2);
  TCanvas* c1 = new TCanvas("c1", "First canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
  c1->Draw();
  c1->Divide(2,2);
  hTrueVtx_X->SetLineColor(kRed);
  hTrueVtx_X->GetXaxis()->SetTitle("True Position X [cm]");
  hTrueVtx_X->GetYaxis()->SetTitle("Num Events");
  hTrueVtx_Y->SetLineColor(kRed);
  hTrueVtx_Y->GetXaxis()->SetTitle("True Position Y [cm]");
  hTrueVtx_Y->GetYaxis()->SetTitle("Num Events");
  hTrueVtx_Z->SetLineColor(kRed);
  hTrueVtx_Z->GetXaxis()->SetTitle("True Position Z [cm]");
  hTrueVtx_Z->GetYaxis()->SetTitle("Num Events");
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
  bsgy->SetLineColor(kBlue);
  bsgy->GetXaxis()->SetTitle("Bonsai Goodness of Fit");
  bsgy->GetYaxis()->SetTitle("Num Events");
  hvtxR->SetLineColor(kBlue);
  hvtxR->GetXaxis()->SetTitle("Reconstructed Event Radius [cm]");
  hvtxR->GetYaxis()->SetTitle("Num Events");
  hVtxDiff_X->GetXaxis()->SetTitle("Reconstructed X Position Error [cm]");
  hVtxDiff_X->GetYaxis()->SetTitle("Num Events");
  hVtxDiff_Y->GetXaxis()->SetTitle("Reconstructed Y Position Error [cm]");
  hVtxDiff_Y->GetYaxis()->SetTitle("Num Events");
  hVtxDiff_Z->GetXaxis()->SetTitle("Reconstructed Z Position Error [cm]");
  hVtxDiff_Z->GetYaxis()->SetTitle("Num Events");
  hNumDigits->GetXaxis()->SetTitle("Num Digits in Event");
  hNumDigits->GetYaxis()->SetTitle("Num Events");
  TLegend *leg1 = new TLegend(0.2033426,0.8192982,0.45961,0.9210526,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextFont(62);
  leg1->SetTextSize(0.04093567);
  leg1->AddEntry(hTrueVtx_X,"True","l");
  leg1->AddEntry(hRecoVtx_X,"Reconstructed","l");
  
  c1->cd(1); hTrueVtx_X->Draw();
  c1->cd(2); hTrueVtx_Y->Draw();
  c1->cd(3); hTrueVtx_Z->Draw();
  
  c1->cd(1); hRecoVtx_X->Draw("same");
  c1->cd(2); hRecoVtx_Y->Draw("same");
  c1->cd(3); hRecoVtx_Z->Draw("same");
  
  c1->cd(1); leg1->Draw();
  
  c1->cd(4); bsgy->Draw();
  
  fileout->cd();
  hTrueVtx_X->Write();
  hTrueVtx_Y->Write();
  hTrueVtx_Z->Write();
  hRecoVtx_X->Write();
  hRecoVtx_Y->Write();
  hRecoVtx_Z->Write();
  hRecoVtx_T->Write();
  bsgy->Write();
  hvtxR->Write();
  hNumDigits->Write();
  hVtxDiff_X->Write();
  hVtxDiff_Y->Write();
  hVtxDiff_Z->Write();
//  fileout->Close();
//  delete fileout;
  return 0;
}
