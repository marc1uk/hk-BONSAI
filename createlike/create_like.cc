#include <stdio.h>
#include <iostream>
#include <math.h>
#include "binfile.h"
#include <stdlib.h>

#define FOURLOG10 9.2103404
#define TBIN 0.4

int main(int argc,char **argv)
{
	float qmin=-1;
	unsigned int nlike=1;
	unsigned int nneg=730;// This value * 0.4 ns is minimum value of likelihood.
	unsigned int bglike=1;

	const unsigned int nbin=840;
	double dx,dhisto[nbin];//={1.0001e-4,6.31e-4,3.981e-3,.02512,.1585,1,.3981,.1585,.0631,.02512,0.01};
	unsigned int histo[nbin];
	binfile      bf(argv[1],'w');
	void         *starts[2];
	int          sizes[3];
	int          numbers[3];
	unsigned int array[3];
	int i;

	double x;
	int j = 0;
	if (argc>2) j = atoi(argv[2]);

	for(i=-((int) nneg); i<((int) nbin)-((int) nneg); i++)
	{
		x=(i*TBIN);
		dhisto[i+nneg]=1.0001e-4;// BG component?

		/*Cylinder 100m*/
		if(j==0){
			std::cout << "/*Cylinder 100m*/" << std::endl;
			dx=(i*TBIN/2.24448);
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);

			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0157*exp(6.76653e-03*x);
				dhisto[i+nneg]+=2.43e-02*exp(4.01711e-02*x);
				dhisto[i+nneg]+=1.37*exp(3.89725e-01*x);
			}
		}

		/*Cylinder 80m*/
		else if (j==1){
			std::cout << "/*Cylinder 80m*/" << std::endl;
			dx=(i*TBIN/2.11462 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0157377*exp(0.00737087*x);
				dhisto[i+nneg]+=1.3722*exp(0.398836*x);
				dhisto[i+nneg]+=0.0233959*exp(0.0436231*x);
			}
		}

		/*Cylinder 60m*/
		else if (j==2){
			std::cout << "/*Cylinder 60m*/" << std::endl;
			dx=(i*TBIN/2.12161 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=1.37383*exp(0.397177*x);
				dhisto[i+nneg]+=0.0164651*exp(0.00752997*x);
				dhisto[i+nneg]+=0.0230223*exp(0.044753*x);
			}
		}

		/*SK 20 B&L*/
		else if (j==3){
			std::cout << "/*SK 20 B&L*/" << std::endl;
			dx=(i*TBIN/1.57723 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0212155*exp(0.090605*x);
				dhisto[i+nneg]+=1.29291*exp(0.449035*x);
				dhisto[i+nneg]+=0.0185637*exp(0.0136284*x);
			}
		}
		/*SK 20 PMT*/
		else if (j==4){
			std::cout << "/*SK 20 PMT*/" << std::endl;
			dx=(i*TBIN/3.62877 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0152849*exp(0.0522208*x);
				dhisto[i+nneg]+=1.4203*exp(0.386938*x);
				dhisto[i+nneg]+=0.0287886*exp(0.0135527*x);
			}
		}
		/*SK 12 B&L*/
		else if (j==5){
			std::cout << "/*SK 20 PMT*/" << std::endl;
			dx=(i*TBIN/1.53266 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0159839*exp(0.0656018*x);
				dhisto[i+nneg]+=1.51575*exp(0.442183*x);
				dhisto[i+nneg]+=0.0162593*exp(0.0133734*x);
			}
		}

		/*HK 20 B&L*/
		else if (j==6){
			dx=(i*TBIN/1.73222 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0271479*exp(0.0623289*x);
				dhisto[i+nneg]+=1.31091*exp(0.429584*x);
				dhisto[i+nneg]+=0.0143397*exp(0.0108354*x);
			}
		}
		/*HK 20 PMT*/
		else if (j==7){
			dx=((i*TBIN)/3.71627 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx);
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0270048*exp(0.0418308*x);
				dhisto[i+nneg]+=1.39011*exp(0.375409*x);
				dhisto[i+nneg]+=0.0206105*exp(0.0102117*x);
			}
		}
		/*Surface detector*/
		else if (j==8){
			dx=(i*TBIN/-1.76557);
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx); 
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=0.0691864*exp(0.0616537*x); 
				double dx2=(x+-8.27888)/14.0146;
				dhisto[i+nneg]+=-0.035261*exp(-0.5*dx2*dx2);
			}
		}
		/* annie */
		else if (j==9){
		// XXX XXX    MAKE SURE THE FORM OF THIS PRINTED IN ANALYSIS_NEW.C   XXX XXX
		// XXX ACTUALLY REPRESENTS THE EXPRESSION AND ALL PARAMETERS OF YOUR FIT XXX
//			// ? which version of wcsim?
//			dx=(i*TBIN/1.77643 );
//			if (x>-1.5) 
//				dhisto[i+nneg]+=1*exp(-0.5*dx*dx); 
//			if (x<=-1.5) 
//			{ 
//				dhisto[i+nneg]+=0.118287*exp(0.165719*x); 
//				dhisto[i+nneg]+=1.32152*exp(1.10514*x); 
//				dhisto[i+nneg]+=3.80723e-48*exp(0.75*x); 
//			}
			// version with 180 LAPPDs and 180 PMTs
			dx=(i*TBIN/-1.7665 );
			if (x>-1.5)
				dhisto[i+nneg]+=1*exp(-0.5*dx*dx); 
			if (x<=-1.5)
			{
				dhisto[i+nneg]+=1.35909*exp(1.08836*x); 
				dhisto[i+nneg]+=0.111633*exp(0.170502*x); 
			}
		}
	}
	sizes[0]=4;
	sizes[1]=4;
	sizes[2]=-1;
	numbers[0]=3;
	numbers[1]=nbin;
	numbers[2]=-1;
	starts[0]=array;
	starts[1]=histo;

	array[0]=nlike;
	array[1]=nneg;
	((float *)starts[0])[2]=qmin;

	for(i=0; i<numbers[1]; i++)
		histo[i]=(unsigned int) (1e5*(1+0.25*log(dhisto[i])/log(10))+0.5);
	//for(i=0; i<numbers[1]; i++)
	//  printf("%u\n",histo[i]);
	bf.write(sizes,numbers,starts);
}

Double_t combinedfunc(Double_t* x, Double_t* par){
	// COPY THIS TEXT FROM ANALYSIS_NEW TO SET PARAMETERS:
	
	// crossover point - upper limit of landaugaus fitrange used in analysis_new.C
	static double crossover = -5;
	// landaugaus parameters printed following fit:
	static std::vector<double> fit_vals{1.,2.,3.,4.};
	
	// printed expression for the exponential part of the fit, including parameters
	
	// END OF COPY SECTION
	
	if((*x)<par[6]) return langaufun(x, par);
	return exp(par[4]+par[5]*(*x));
}

