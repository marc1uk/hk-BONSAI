void analysis(){
	//TFile *file0 = TFile::Open("Cylinder100_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("Cylinder80_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("Cylinder60_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("SK_20BL_10MeV_Uni_Iso.tune.root");
	TFile *file0 = TFile::Open("SK_20PMT_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("SK_12BL_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("HK20BL_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("HK20PMT_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("Cylinder8x14_10MeV_Uni_Iso.tune.root");
	//TFile *file0 = TFile::Open("SKtest_10MeV_Uni_Iso.tune.root");

	TF1 *fit = new  TF1();
	//TF1 *fit = new TF1("fit","gaus(0)*(x<1.5)+(x>=1.5)*(expo(3)+gaus(5)+expo(8))",-60,250);
	TF1 *fit = new TF1("fit","gaus(0)*(x<1.5)+(x>=1.5)*(expo(3)+expo(5)+expo(7))",-60,250);
	fit->SetNpx(1000);
	//expo = exp([0]+[1]*x)
	likelihood->Draw();
	fit ->SetParameter(0,5.85655e+04);
	fit ->SetParameter(1,4.72705e-01);
	fit ->SetParameter(2,2.68512e+00);
	fit ->SetParameter(3,7.53496e+00);
	fit ->SetParameter(4,-9.49172e-03);
	fit ->SetParameter(5,1.05580e+01);
	fit ->SetParameter(6,-2.81987e-01);
	/*
   1  p0           6.56779e+04   9.50418e+01  -3.36903e+00   6.80090e-03
   2  p1          -5.42404e-02   6.01265e-03   1.11805e-02  -1.50516e+02
   3  p2           2.12371e+00   4.19033e-03   6.25285e-03   4.88323e+02
   4  p3           7.30748e+00   3.94994e-03   1.70887e-03   3.82046e+03
   5  p4          -9.49172e-03   3.57330e-05   4.09642e-05   3.41234e+05
   6  p5           1.05580e+01   2.26308e-02  -9.71021e-02   4.04517e+03
   7  p6          -2.81987e-01   2.28885e-03   8.15921e-03   1.76097e+04
   8  p7           1.10558e+01   2.63455e-02   2.49574e-02   1.81210e+03
   9  p8          -5.62798e-01   6.30956e-03   3.22930e-02   5.46450e+03
	   fit ->SetParameter(5,4.14563e+07);
	   fit ->SetParameter(6,-4.05867e+01);
	   fit ->SetParameter(7,1.14106e+01);
	   */
	likelihood->Fit("fit");
	  fit ->SetParameter(7,0);
	      fit ->SetParameter(8,0);
	likelihood->Fit("fit");

	  fit ->SetParameter(7,0);
	      fit ->SetParameter(8,0);
	likelihood->Fit("fit");

	//cout <<" dx=((i*TBIN+" << fit->GetParameter(1) << ")/" << fit->GetParameter(2) << " );" << endl;
	cout <<" dx=(i*TBIN/" << fit->GetParameter(2) << " );" << endl;
	cout <<" if (x>-1.5) "<< endl;
	cout <<" dhisto[i+nneg]+=1*exp(-0.5*dx*dx); "<< endl;
	cout <<" if (x<=-1.5) "<< endl;
	cout <<" { "<< endl;
	cout <<" dhisto[i+nneg]+=" << exp(fit->GetParameter(3))/fit->GetParameter(0) << "*exp(" << -1*(fit->GetParameter(4)) << "*x); "<< endl;
	cout <<" dhisto[i+nneg]+=" << exp(fit->GetParameter(5))/fit->GetParameter(0) << "*exp(" << -1*(fit->GetParameter(6)) << "*x); "<< endl;
	cout <<" dhisto[i+nneg]+=" << exp(fit->GetParameter(7))/fit->GetParameter(0) << "*exp(" << -1*(fit->GetParameter(8)) << "*x); "<< endl;
	cout <<" } "<< endl;
}
