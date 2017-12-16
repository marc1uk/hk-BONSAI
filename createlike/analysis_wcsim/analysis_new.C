#include <iostream>
#include <TSystem.h>
#include <TH1F.h>
#include <TF1.h>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TStyle.h>
#include <TFile.h>
#include <math.h>
#include <TMath.h>
#include "../../gitver/Bonsai_v0/bonsai/pmt_geometry.h"
#include "../../gitver/Bonsai_v0/bonsai/likelihood.h"
#include "../../gitver/Bonsai_v0/bonsai/goodness.h"

#include "WCSimRootGeom.hh" // << ENSURE THE CORRECT VERSIONS ARE CHECKED OUT
#include "WCSimRootEvent.hh"
//#include "WCSimRootLinkDef.hh"

#include <regex>

// defined in $ROOTSYS/tutorials/fit/langaus.C
//#include "langaus.C"
// function that fits the parameters
TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
// function which represents the actual mapping of x->y (the 'fit formula')
Double_t langaufun(Double_t *x, Double_t *par);
Double_t combinedfunc(Double_t* x, Double_t* par); // defined at the bottom of this file

namespace{
    constexpr double triggeroffset = 0;//950.0; // WCSimRootTrigger::offset. 
    // TODO: put this into WCSimRootOptions so it isn't hard-coded!
    constexpr bool ZERO_DIGIT_TS = true; // remove triggeroffset and trigger time to zero digit times
    
    constexpr double LIGHT_SPEED = 21.58333; // units?
    
    // CHANGE IF ANALYZING OTHER DETECTORS
    constexpr bool isANNIE = true; // Swap Y and Z axes, along with any other required changes.
}

void analysis_new(char *filename=NULL, bool surfaceDetector=false, bool verbose=false){
//    char* wcsimdirenv;
//    wcsimdirenv = getenv ("WCSIMDIR");
//    if(wcsimdirenv !=  NULL){
//          gSystem->Load("${WCSIMDIR}/wcsim/libWCSimRoot.so");
//          gSystem->SetIncludePath("${WCSIMDIR}/wcsim/include");
//    }else{
//        gSystem->Load("../libWCSimRoot.so");
//    }
    gROOT->Reset();
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    gSystem->Load("${WCSIMDIR}/wcsim/libWCSimRoot.so");
    gSystem->SetIncludePath("${WCSIMDIR}/wcsim/include");

    TFile *file;
    // Open the file
    if (filename==NULL){
        file = new TFile("../wcsim.root","read");
    }else{
        file = new TFile(filename,"read");
    }
    if (!file->IsOpen()){
        std::cout << "Error, could not open input file: " << filename << std::endl;
        return;
    }

    TString infile = filename;

    TH1F *likelihood;

    if(file->GetListOfKeys()->Contains("wcsimT")) {

        TString outfile = infile.Insert(infile.Length()-5,".tune");
        TFile *file_out = new TFile(outfile,"recreate");

        // Get the a pointer to the tree from the file
        TTree *tree = (TTree *) file->Get("wcsimT");

        // Get the number of events
        int nevent = tree->GetEntries();
        if (verbose) printf("nevent %d\n", nevent);

        // Create a WCSimRootEvent to put stuff from the tree in

        WCSimRootEvent *wcsimrootsuperevent = new WCSimRootEvent();

        // Set the branch address for reading from the tree
        TBranch *branch = tree->GetBranch("wcsimrootevent");
        branch->SetAddress(&wcsimrootsuperevent);

        // Force deletion to prevent memory leak
        tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

        // Geometry tree - only need 1 "event"
        TTree *geotree = (TTree *) file->Get("wcsimGeoT");
        WCSimRootGeom *geo = 0;
        geotree->SetBranchAddress("wcsimrootgeom", &geo);
        if (verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
        if (geotree->GetEntries() == 0) {
            exit(9);
        }
        geotree->GetEntry(0);

        // start with the main "subevent", as it contains most of the info
        // and always exists.
        WCSimRootTrigger *wcsimrootevent;

        TH1F *h1 = new TH1F("PMT Hits", "PMT Hits", 500, 0, 500);
        TH1F *hvtx01 = new TH1F("Event VTX X", "Event VTX X", 200, -1500, 1500);
        TH1F *hvtx02 = new TH1F("Event VTX Y", "Event VTX Y", 200, -1500, 1500);
        TH1F *hvtx03 = new TH1F("Event VTX Z", "Event VTX Z", 200, -1500, 1500);
        likelihood = new TH1F("likelihood", "t-tof", 750, -25, 50);

        // Now loop over events
        for (int ev = 0; ev < nevent; ev++) {
            // Read the event from the tree into the WCSimRootEvent instance
            tree->GetEntry(ev);
            wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
            //if(verbose){
            float vtxX, vtxY, vtxZ;
            if(not isANNIE){
              vtxX = wcsimrootevent->GetVtx(0)-geo->GetWCOffset(0); // offsets added 13-12-17
              vtxY = wcsimrootevent->GetVtx(1)-geo->GetWCOffset(1);
              vtxZ = wcsimrootevent->GetVtx(2)-geo->GetWCOffset(2);
            } else {
              vtxX = wcsimrootevent->GetVtx(0)-geo->GetWCOffset(0); // offsets added 13-12-17;
              vtxY = wcsimrootevent->GetVtx(2)-geo->GetWCOffset(2);;
              vtxZ = wcsimrootevent->GetVtx(1)-geo->GetWCOffset(1);;
            }
            cout<<"event "<<ev<<" had e- vertex at ("<<vtxX<<","<<vtxY<<","<<vtxZ<<")"<<endl;
            hvtx01->Fill(vtxX);
            hvtx02->Fill(vtxY);
            hvtx03->Fill(vtxZ);

            int i;

            // Now look at the Cherenkov hits

            // Get the number of Cherenkov hits.
            // Note... this is *NOT* the number of photons that hit tubes.
            // It is the number of tubes hit with Cherenkov photons.
            // The number of digitized tubes will be smaller because of the threshold.
            // Each hit "raw" tube has several photon hits.  The times are recorded.
            // See http://nwg.phy.bnl.gov/DDRD/cgi-bin/private/ShowDocument?docid=245
            // for more information on the structure of the root file.
            //

            int ncherenkovhits = wcsimrootevent->GetNcherenkovhits();
            cout<<"it had "<<ncherenkovhits<<" cherenkov hits"<<endl;

            // Grab the big arrays of times and parent IDs
            TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();

            // Look at digitized hit info

            // Get the number of digitized hits
            // Loop over sub events
            if (verbose) std::cout << "DIGITIZED HITS:" << std::endl;
            for (int index = 0; index < wcsimrootsuperevent->GetNumberOfEvents(); index++) {
                wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
                WCSimRootEventHeader* header = wcsimrootevent->GetHeader();
                int triggertime = header->GetDate();
                if (verbose) std::cout << "Sub event number = " << index << "\n";

                int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
                if (ncherenkovdigihits > 0) h1->Fill(ncherenkovdigihits);
                if (verbose){ printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
                cout<<" and "<<ncherenkovdigihits<<" digi hits"<<endl; }

                for (i = 0; i < ncherenkovdigihits; i++) {
                    // Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
                    TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
                    WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit *>(element);
                    float T = wcsimrootcherenkovdigihit->GetT();
                    if(ZERO_DIGIT_TS) T += triggertime - triggeroffset;
                    
                    int pmtID = wcsimrootcherenkovdigihit->GetTubeId();
                    WCSimRootPMT pmt = geo->GetPMT(pmtID - 1);
                    float pmtX, pmtY, pmtZ;
                    if(not isANNIE){
                      pmtX = pmt.GetPosition(0);
                      pmtY = pmt.GetPosition(1);
                      pmtZ = pmt.GetPosition(2);
                    } else {
                      pmtX = pmt.GetPosition(0);
                      pmtY = pmt.GetPosition(2);
                      pmtZ = pmt.GetPosition(1);
                    }
                    float ttof = T - TMath::Sqrt((vtxX - pmtX) * (vtxX - pmtX) + (vtxY - pmtY) * (vtxY - pmtY) +
                                                 (vtxZ - pmtZ) * (vtxZ - pmtZ)) / LIGHT_SPEED;
                    likelihood->Fill(ttof);
                } // End of loop over Cherenkov digihits
            } // End of loop over trigger
            // reinitialize super event between loops.
            wcsimrootsuperevent->ReInitialize();
        } // End of loop over events
        file_out->Write();
    }
    else if(file->GetListOfKeys()->Contains("likelihood")){
        likelihood = (TH1F*) file->Get("likelihood");
    }
    else{
        cout << "Error: Input file is not wcsim output and does not contain likelihood histogram.";
        return;
    }

    if(!surfaceDetector) {	//use this for annie
//        TF1 *fit = new TF1("fit", "gaus(0)*(x<1.5)+(x>=1.5)*(expo(3)+expo(5)+expo(7))", -60, 250);
//        TF1 *fit = new TF1("fit", "gaus(0)*(x<1.5)+(x>=1.5)*(expo(3)+expo(5))", -60, 250); // 180 LAPPD version
//        TF1 *fit = new TF1("fit", "gaus(0)*(x<[3])+(x>[4])*(expo(5))+gaus(7)+((x>[10])*[11])*expo(12)"); // 180 PMTs of varying size version .... has issues
        
        ///////////////////
        // Landau gaus fit... works well, but is much more awkward. :(
         // fit range, lower limit, upper limit
        std::vector<double> fitrange{-25.,-0.3};  // UPPERLIMIT HERE IS THE INTERFACE PARAMETER XXX XXX XXX XXX XXX
        // starting fit vals: landau width, landau centre, integral, width of gaussian component
        std::vector<double> start_vals{2.,-6.,likelihood->Integral(),1.6};
        // lower limits on fit parameter ranges
        std::vector<double> param_lowerlims{0.,-100.,1.,0.01};
        // upper limits on fit parameter ranges
        std::vector<double> param_upperlims{10.,100.,100.*likelihood->Integral(),100.};
        // fit results
        std::vector<double> fit_vals(4), fit_val_errors(4);
        Double_t chisqr;                  // fit chi^2
        Int_t    ndf;                     // fit NDF

        // first fit the former half of the curve with a landau function convoluted with gaussian
        // this provides better matching at the interface between the exponential latter half
        // and the gaussian front.
        TF1 *fitlangaus = langaufit(likelihood,fitrange.data(),start_vals.data(),param_lowerlims.data(),
                              param_upperlims.data(),fit_vals.data(),fit_val_errors.data(),&chisqr,&ndf);
        
        //TF1* fitexp = new TF1("fitexp",TString::Format("(x>%f)*expo(0)",fitrange.at(1)),fitrange.at(1),50.);
        //TF1* fit = new TF1("fit", "Fitfcn_likelihood+fitexp"); // first one is name of langaus fn.
        // ^ can't combine like this when one TF1 is a custom TFormulae/c++ function. ^
        // instead we need to define another custom c++ function. 
        
        // fit the second half, an exponential.
        // the only reasonable way to do this is to hard-code it in combinedfunc.
        // the function is just 'expo(0)' equivalent to 'exp(par[0]+par[1]*x)'
        std::vector<double> exp_vals{6.5,-0.2}; // initial fit parameters
        
        // combine the fit parameters from the two parts.
        std::vector<double> allparams(fit_vals);
        allparams.insert(allparams.end(),exp_vals.begin(),exp_vals.end());
        allparams.push_back(fitrange.at(1));
        
        TF1* fit = new TF1("fit",combinedfunc,fitrange.at(0),50.,allparams.size());
        fit->SetParameters(allparams.data());
        fit->SetNpx(1000);
        likelihood->Fit("fit","MR");
        
        TCanvas* c = (TCanvas*)gROOT->FindObject("c1");
        c->SetLogy();
        c->Clear();
        likelihood->Draw();
        likelihood->Draw("lsame");
        c->Update();
        
//        a good set of fit parameters:
//           1  p0           2.38468e+00   4.16930e-02   1.02382e-04   4.19183e-06
//           2  p1          -6.55417e+00   2.91117e-02   1.47047e-04  -3.70137e-06
//           3  p2           2.02002e+04   1.68480e+02   7.17281e-01  -8.17844e-10
//           4  p3           1.60740e+00   4.68113e-02   1.85193e-04  -1.11266e-05
//           5  p4           6.58271e+00   7.26704e-03   6.36750e-05   1.52431e-05
//           6  p5          -2.06821e-01   1.14914e-03   1.01137e-05  -1.05242e-04
//           7  p6          -4.25588e-01   1.17709e-01   1.62538e-01   1.14312e-01

        
        ///////////////////
//        TCanvas * c = new TCanvas();
//        c->SetLogy();
//        likelihood->Draw();
//        TF1 *fit = new TF1("fit", "gaus(0)*(x<[3])+(x>[3])*(expo(4))+gaus(6)+gaus(13)");
//        fit->SetNpx(1000);
////        // bonsai default?
//////        fit->SetParameter(0, 5.85655e+04);
////        fit->SetParameter(1, 4.72705e-01);
////        fit->SetParameter(2, 2.68512e+00);
////        fit->SetParameter(3, 7.53496e+00);
////        fit->SetParameter(4, -9.49172e-03);
////        fit->SetParameter(5, 1.05580e+01);
////        fit->SetParameter(6, -2.81987e-01);
////        fit->SetParameter(1, 1000.);
////        // 180 LAPPDs+PMTs
////        fit->SetParameter(2, -2.);
////        fit->SetParameter(3, 2.);
////        fit->SetParameter(4, -9.49172e-03);
////        fit->SetParameter(5, 1.05580e+01);
////        fit->SetParameter(6, -2.81987e-01);
////        likelihood->Fit("fit");
////        fit->SetParameter(7, -100);	//modified from 0 - annie
////        fit->SetParameter(8, 0);
////        likelihood->Fit("fit");
////        fit->SetParameter(7, -100);	// modified from 0 - annie
////        fit->SetParameter(8, 0);
//        // 124 PMTs only, varying sizes
//        fit->SetParameter(1, -6.6);
//        fit->SetParameter(2, -2.7);
//        fit->SetParameter(3, -5.);
//        fit->SetParameter(4, 6.4);
//        fit->SetParameter(5, -0.18);
//        fit->SetParameter(6, 300.);
//        fit->SetParameter(7, -3.);
//        fit->SetParameter(8, 5.);
//        fit->SetParameter(9, -7.);
//        fit->SetParameter(10, 1000.);
//        fit->SetParameter(11, 0.5);
//        fit->SetParameter(12, -5.e-03);
//        fit->SetParameter(13, 4.5);
//        fit->SetParameter(14, 14.);
//        fit->SetParameter(15, 4.);
//        likelihood->Fit("fit");
//        cout << "  dx=(i*TBIN/" << fit->GetParameter(2) << " );" << endl;
//        cout << "  if (x>-1.5) " << endl;
//        cout << "  dhisto[i+nneg]+=1*exp(-0.5*dx*dx); " << endl;
//        cout << "  if (x<=-1.5) " << endl;
//        cout << "  { " << endl;
//        cout << "    dhisto[i+nneg]+=" << exp(fit->GetParameter(3)) / fit->GetParameter(0) << "*exp(" << -1 * (fit->GetParameter(4)) << "*x); " << endl;
//        cout << "    dhisto[i+nneg]+=" << exp(fit->GetParameter(5)) / fit->GetParameter(0) << "*exp(" << -1 * (fit->GetParameter(6)) << "*x); " << endl;
//        //cout << "    dhisto[i+nneg]+=" << exp(fit->GetParameter(7)) / fit->GetParameter(0) << "*exp(" << -1 * (fit->GetParameter(8)) << "*x); " << endl;
//        cout << "  } " << endl;
    }
    else {
        TF1 *fit = new TF1("fit", "gaus(0)+(x>=2)*(expo(3)+gaus(5))", -60, 250);
        // guas: aplitude, mean, width
        fit->SetNpx(1000);
        likelihood->Draw();
//        fit->SetParameter(0, 5.85655e+04);
//        fit->SetParameter(1, 4.72705e-01);
//        fit->SetParameter(2, 2.68512e+00);
//        fit->SetParameter(3, 7.53496e+00);
//        fit->SetParameter(4, -9.49172e-03);
//        fit->SetParameter(5, 5.85655e+06);
//        fit->SetParameter(6, 10);
//        fit->SetParameter(7, 10);
        fit->SetParameter(0, 1000.);
        fit->SetParameter(1, -1.);
        fit->SetParameter(2, 2.);
        fit->SetParameter(3, 5.);
        fit->SetParameter(4, -0.1);
        fit->SetParameter(5, 1000.);
        fit->SetParameter(6, 10.);
        fit->SetParameter(7, 10.);
        likelihood->Fit("fit");
        cout << "Test4" << endl;
        likelihood->Fit("fit");
        likelihood->Fit("fit");
        cout << "  dx=(i*TBIN/" << fit->GetParameter(2) << " );" << endl;
        cout << "  dhisto[i+nneg]+=1*exp(-0.5*dx*dx); " << endl;
        cout << "  if (x<=-1.5) " << endl;
        cout << "  { " << endl;
        cout << "    dhisto[i+nneg]+=" << exp(fit->GetParameter(3)) / fit->GetParameter(0) << "*exp(" << -1 * (fit->GetParameter(4)) << "*x); " << endl;
        cout << "    double dx2=(x+" << fit->GetParameter(6) << ")/" << fit->GetParameter(7) << ";" << endl;
        cout << "    dhisto[i+nneg]+=" << fit->GetParameter(5) / fit->GetParameter(0) << "*exp(-0.5*dx2*dx2);" << endl;
        cout << "  } " << endl;
    }

}

Double_t combinedfunc(Double_t* x, Double_t* par){
    if((*x)<par[6]) return langaufun(x, par);
    return exp(par[4]+par[5]*(*x));
}
