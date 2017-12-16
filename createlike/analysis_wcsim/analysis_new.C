#include <iostream>
#include <TH1F.h>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TTree.h>
#include <TCanvas.h>
//#include <../include/WCSimRootEvent.hh>
//#include <../include/WCSimRootGeom.hh>
//#include <../include/WCSimRootLinkDef.hh>

void analysis_new(char *filename=NULL, bool surfaceDetector=false, bool verbose=false){
    char* wcsimdirenv;
    wcsimdirenv = getenv ("WCSIMDIR");
    if(wcsimdirenv !=  NULL){
        gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
    }else{
        gSystem->Load("../libWCSimRoot.so");
    }

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
        likelihood = new TH1F("likelihood", "t-tof", 750, -50, 250);

        // Now loop over events
        for (int ev = 0; ev < nevent; ev++) {
            // Read the event from the tree into the WCSimRootEvent instance
            tree->GetEntry(ev);
            wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
            //if(verbose){
            float vtxX = wcsimrootevent->GetVtx(0);
            float vtxY = wcsimrootevent->GetVtx(1);
            float vtxZ = wcsimrootevent->GetVtx(2);
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

            // Grab the big arrays of times and parent IDs
            TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();

            // Look at digitized hit info

            // Get the number of digitized hits
            // Loop over sub events
            if (verbose) std::cout << "DIGITIZED HITS:" << std::endl;
            for (int index = 0; index < wcsimrootsuperevent->GetNumberOfEvents(); index++) {
                wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
                if (verbose) std::cout << "Sub event number = " << index << "\n";

                int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
                if (ncherenkovdigihits > 0) h1->Fill(ncherenkovdigihits);
                if (verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);

                for (i = 0; i < ncherenkovdigihits; i++) {
                    // Loop through elements in the TClonesArray of WCSimRootCherenkovDigHits
                    TObject *element = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
                    WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit *>(element);
                    float T = wcsimrootcherenkovdigihit->GetT();
                    int pmtID = wcsimrootcherenkovdigihit->GetTubeId();
                    WCSimRootPMT pmt = geo->GetPMT(pmtID - 1);
                    float pmtX = pmt.GetPosition(0);
                    float pmtY = pmt.GetPosition(1);
                    float pmtZ = pmt.GetPosition(2);
                    float ttof = T - TMath::Sqrt((vtxX - pmtX) * (vtxX - pmtX) + (vtxY - pmtY) * (vtxY - pmtY) +
                                                 (vtxZ - pmtZ) * (vtxZ - pmtZ)) / 21.58333;
                    likelihood->Fill(ttof);
                } // End of loop over Cherenkov digihits
            } // End of loop over trigger
            // reinitialize super event between loops.
            wcsimrootsuperevent->ReInitialize();
        } // End of loop over events
        file_out->Write();
    }
    else if(file->GetListOfKeys()->Contains("likelihood")){
        likelihood = (TH1F*) infile->Get("likelihood");
    }
    else{
        cout << "Error: Input file is not wcsim output and does not contain likelihood histogram.";
        return -1;
    }

    if(!surfaceDetector) {
        TF1 *fit = new TF1("fit", "gaus(0)*(x<1.5)+(x>=1.5)*(expo(3)+expo(5)+expo(7))", -60, 250);
        fit->SetNpx(1000);
        TCanvas * c = new TCanvas();
        c->SetLogy();
        likelihood->Draw();
        fit->SetParameter(0, 5.85655e+04);
        fit->SetParameter(1, 4.72705e-01);
        fit->SetParameter(2, 2.68512e+00);
        fit->SetParameter(3, 7.53496e+00);
        fit->SetParameter(4, -9.49172e-03);
        fit->SetParameter(5, 1.05580e+01);
        fit->SetParameter(6, -2.81987e-01);
        likelihood->Fit("fit");
        fit->SetParameter(7, 0);
        fit->SetParameter(8, 0);
        likelihood->Fit("fit");
        fit->SetParameter(7, 0);
        fit->SetParameter(8, 0);
        likelihood->Fit("fit");
        cout << "  dx=(i*TBIN/" << fit->GetParameter(2) << " );" << endl;
        cout << "  if (x>-1.5) " << endl;
        cout << "  dhisto[i+nneg]+=1*exp(-0.5*dx*dx); " << endl;
        cout << "  if (x<=-1.5) " << endl;
        cout << "  { " << endl;
        cout << "    dhisto[i+nneg]+=" << exp(fit->GetParameter(3)) / fit->GetParameter(0) << "*exp(" << -1 * (fit->GetParameter(4)) << "*x); " << endl;
        cout << "    dhisto[i+nneg]+=" << exp(fit->GetParameter(5)) / fit->GetParameter(0) << "*exp(" << -1 * (fit->GetParameter(6)) << "*x); " << endl;
        cout << "    dhisto[i+nneg]+=" << exp(fit->GetParameter(7)) / fit->GetParameter(0) << "*exp(" << -1 * (fit->GetParameter(8)) << "*x); " << endl;
        cout << "  } " << endl;
    }
    else {
        TF1 *fit = new TF1("fit", "gaus(0)+(x>=2)*(expo(3)+gaus(5))", -60, 250);
        fit->SetNpx(1000);
        likelihood->Draw();
        fit->SetParameter(0, 5.85655e+04);
        fit->SetParameter(1, 4.72705e-01);
        fit->SetParameter(2, 2.68512e+00);
        fit->SetParameter(3, 7.53496e+00);
        fit->SetParameter(4, -9.49172e-03);
        fit->SetParameter(5, 5.85655e+06);
        fit->SetParameter(6, 10);
        fit->SetParameter(7, 10);
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
