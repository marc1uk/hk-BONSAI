// Based on Root test Event.cxx
////////////////////////////////////////////////////////////////////////
#ifndef WCSIMBONSAIVERBOSE
//#define WCSIMBONSAIVERBOSE 1
#endif

//#include "G4ios.hh"
#include "TObject.h"
#include "TDirectory.h"
#include "WCSimBonsai.hh"

#include <stdio.h>
#include "vertexfit.h"
#include "pmt_geometry.h"
#include "fourhitgrid.h"
#include "goodness.h"
#include "likelihood.h"
#include "binfile.h"
//#include "bscalls.h"
#include <TMath.h>
#include <iostream>

vertexfit *vf_root;
likelihood *tf_root;
pmt_geometry *geom_root;
binfile    *bsfile_root;
comtype pmt_array;

#define MINHIT 4

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimBonsai)
#endif

using namespace std;

	//______________________________________________________________________________
WCSimBonsai::WCSimBonsai()
{}

//______________________________________________________________________________
WCSimBonsai::~WCSimBonsai()
{
	if(tf_root!=NULL)tf_root->set_hits(NULL);
	delete(tf_root);
	delete(geom_root);
}

Int_t WCSimBonsai::Init(WCSimRootGeom *fGeo, bool swapYandZ=false){

//TODO use fOrientation from fGeo to avoid the need for this. Right now it's not properly set in WCSim though
// n.b it's set to 0 in the RunAction ...?
// other consideraitons: use fgeo_type (an int) to set detector-specific configurations.

		int fNPMT = fGeo->GetWCNumPMT();
		// fill the arrays of geometry info
		WCSimRootPMT pmt;

#ifdef WCSIMBONSAIVERBOSE
		std::cout << " Geometry Loading (NPMT=" << fNPMT << ") " << std::endl;
#endif
		double max_cylinder_height=0, max_cylinder_radius=0; //PMT geometory array for bonsai. T. Yano
		for (int ipmt=0; ipmt<fNPMT; ipmt++){    
			pmt = fGeo->GetPMT(ipmt);
			//std::cout<<"<WCSimBonsai::ReadGeom> ipmt="<<ipmt<<" tubeno="<<pmt.GetTubeNo()<<std::endl;

//PMT geometory array for bonsai. T. Yano
			if(swapYandZ){
			// TANK ORIENTATION IS DIFFERENT SO SWAP AROUND INDICES.
			// TANK IS ALSO OFFSET FROM (0,0,0): tank -152<X<152, -212<Y<183, 15<Z<320 cm
			pmt_array.pmt_position[ipmt][0] = pmt.GetPosition(0)-fGeo->GetWCOffset(0); // x
			pmt_array.pmt_position[ipmt][2] = pmt.GetPosition(1)-fGeo->GetWCOffset(1); // y
			pmt_array.pmt_position[ipmt][1] = pmt.GetPosition(2)-fGeo->GetWCOffset(2); // z
			pmt_array.pmt_direction[ipmt][0] = pmt.GetOrientation(0);
			pmt_array.pmt_direction[ipmt][2] = pmt.GetOrientation(1);
			pmt_array.pmt_direction[ipmt][1] = pmt.GetOrientation(2);
			pmt_array.pmt_r[ipmt] = 
			TMath::Sqrt( TMath::Power(pmt.GetPosition(0)-fGeo->GetWCOffset(0),2) + TMath::Power(pmt.GetPosition(2)-fGeo->GetWCOffset(2),2) );
			} else {
			pmt_array.pmt_position[ipmt][0] = pmt.GetPosition(0)-fGeo->GetWCOffset(0);
			pmt_array.pmt_position[ipmt][1] = pmt.GetPosition(1)-fGeo->GetWCOffset(1);
			pmt_array.pmt_position[ipmt][2] = pmt.GetPosition(2)-fGeo->GetWCOffset(2);
			pmt_array.pmt_direction[ipmt][0] = pmt.GetOrientation(0);
			pmt_array.pmt_direction[ipmt][1] = pmt.GetOrientation(1);
			pmt_array.pmt_direction[ipmt][2] = pmt.GetOrientation(2);
			pmt_array.pmt_r[ipmt] = TMath::Sqrt( pmt.GetPosition(0)*pmt.GetPosition(0) + (pmt.GetPosition(1))*(pmt.GetPosition(1)) );
			}
			pmt_array.pmt_type[ipmt] = 1;
			pmt_array.pmt_theta[ipmt] = 0;
			if (max_cylinder_height < pmt_array.pmt_position[ipmt][2]) max_cylinder_height = pmt_array.pmt_position[ipmt][2];
			if (max_cylinder_radius < pmt_array.pmt_r[ipmt]) max_cylinder_radius = pmt_array.pmt_r[ipmt];
			//PMT geometory array for bonsai. T. Yano

#ifdef WCSIMBONSAIVERBOSE
			std::cout << pmt_array.pmt_position[ipmt][0] << " " << pmt_array.pmt_position[ipmt][1] << " " << pmt_array.pmt_position[ipmt][1] << " (" << ipmt/(double)fNPMT << "%)"<< std::endl;
#endif
		}
		//PMT geometory array for bonsai. T. Yano
		pmt_array.n_pmts=fNPMT;
		geom_root = new pmt_geometry(&pmt_array);
		tf_root=new likelihood(max_cylinder_radius, max_cylinder_height);
		return 0;
}

//______________________________________________________________________________
Int_t WCSimBonsai::BonsaiFit(float *vert,float *result,float *maxlike,int *nsel,
		int *nhit,int *cab,float *t,float *q)
{
	// vert[4]   : vertex (x,y,z,t)
	// result[5] : direction (x,y,z,t,ll0) 
	// maxlike[3]: cf.maxq(), tgood, goodness
	// nsel      : bshits.nselected()
	goodness    bshits(tf_root->sets(),tf_root->chargebins(),geom_root,*nhit,cab,t,q);
	*nsel=bshits.nselected();
	if (bshits.nselected()<MINHIT) return(0);
	fourhitgrid gr(geom_root->cylinder_radius(),geom_root->cylinder_height(),&bshits);
	bonsaifit   cf(tf_root);
	int         nfit;
	float       gdn[tf_root->sets()];

	tf_root->set_hits(&bshits);
	tf_root->maximize(&cf,&gr);
	if ((nfit=tf_root->nfit())==0) return(0);
	*vert=cf.xfit();
	vert[1]=cf.yfit();
	vert[2]=cf.zfit();
	maxlike[2]=tf_root->goodness(*maxlike,vert,gdn);
	tf_root->tgood(vert,0,maxlike[1]);
	nsel[1]=tf_root->nwind(vert,-6,12);
	*maxlike=cf.maxq();
	cf.fitresult();
	vert[3]=tf_root->get_zero();
	tf_root->get_dir(result);
	result[5]=tf_root->get_ll0();
	tf_root->set_hits(NULL);
	return(nfit);
}
