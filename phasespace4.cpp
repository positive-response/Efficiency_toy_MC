#include "TMath.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <iostream>

double getGeometricalEfficieny(const double theta_min = 1.0472, const double theta_max = 2.0944, const int minimumNumberOfGamma = 3)
{
const  double electron_mass = 0.000511; // in GeV                                                                                                 
TLorentzVector electron(0.0,0.0,0.0,electron_mass);
TLorentzVector proton(0.0,0.0,0.0,electron_mass);
TLorentzVector parent_particle = electron + proton;
const int numberOfDaughterParticles = 4;
double masses[numberOfDaughterParticles] = {0.0};

TGenPhaseSpace event;
event.SetDecay(parent_particle, numberOfDaughterParticles, masses);

int iteration = 1000000;
double weight = 0.0;
double gamma_theta[numberOfDaughterParticles]={0.0};
int total = 0;
int accepted = 0;
int n_min = 0;
for(int i = 0; i < iteration; i++)
   {
	   n_min = 0;
	   total++;
	   weight = event.Generate();

	   for(int j = 0; j< numberOfDaughterParticles; j++)
	   {
		   TLorentzVector *gamma = event.GetDecay(j);
                   gamma_theta[j] = gamma->Theta();
		  if ((gamma_theta[j] > theta_min) & (gamma_theta[j] < theta_max))
	     		  n_min++;
	   }
	     if (n_min >= minimumNumberOfGamma)
		     accepted++;
	     n_min = 0;
  }

//std::cout<<"ratio: "<< static_cast<double>(accepted)/total<<std::endl;
//std::cout<<"total: "<< total<<std::endl;
//std::cout<<"accepted: "<< accepted<<std::endl;
return static_cast<double>(accepted)/total;
}

int main()
{
	const double mininumTheta = 1.0472;
	const double maximumTheta = 2.0944;
	const int minimumNumberOfGamma = 2;
  std::cout << getGeometricalEfficieny(mininumTheta, maximumTheta, minimumNumberOfGamma) << std::endl;;
}
