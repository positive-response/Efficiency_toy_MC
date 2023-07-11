#include "TMath.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <iostream>

void phasespace4()
{

const  double melectron = 0.000511; // in GeV                                                                                                 
TLorentzVector e(0.0,0.0,0.0,melectron);
TLorentzVector p(0.0,0.0,0.0,melectron);
TLorentzVector W = e + p;
const int npar = 4;
double masses[npar] = {0.0};

TGenPhaseSpace event;
event.SetDecay(W,npar,masses);

int iter = 1000000;
double weight = 0.0;
double gamma_E[npar]={0.0};
double gamma_theta[npar]={0.0};
int total = 0;
int accepted = 0;
int n_min = 0;

for(int i = 0; i < iter; i++)
   {
	   n_min = 0;
	   total++;
	   weight = event.Generate();

	   for(int j = 0; j< npar; j++)
	   {
		   TLorentzVector *gamma = event.GetDecay(j);
                   gamma_theta[j] = gamma->Theta();
		  if (gamma_theta[j] > 1.0472 & gamma_theta[j] < 2.0944)
	     		  n_min++;
	   }
	     if (n_min >= 3)
		     accepted++;
	     n_min = 0;
  }

std::cout<<"ratio: "<< static_cast<double>(accepted)/total<<std::endl;
std::cout<<"total: "<< total<<std::endl;
std::cout<<"accepted: "<< accepted<<std::endl;

}
