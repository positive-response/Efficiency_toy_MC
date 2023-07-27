#include "TMath.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <iostream>

using namespace std;

double getGeometricalEfficieny(const double , const double, const int);
double* calculateAcceptance(const double, const double);

int main()
{
	const double radiusOfDetector = 42.5;  // in cm 
	const double half_lengthOfDetector = 25.0;  //in cm
	double * thetaLimits = calculateAcceptance(radiusOfDetector, half_lengthOfDetector);
	const double mininumTheta = *(thetaLimits + 0);
	const double maximumTheta = *(thetaLimits + 1);
	const int minimumNumberOfGamma = 3;
        std::cout << getGeometricalEfficieny(mininumTheta, maximumTheta, minimumNumberOfGamma) << std::endl;
}

/************************************************************************************************************/
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
return static_cast<double>(accepted)/total;
}
/*****************************************************************************************/
double* calculateAcceptance(const double radius = 42.5, const double half_length = 25.0)
{

	const double minimumAcceptance = TMath::ATan(radius/half_length);
	const double maximumAcceptance = 2 * minimumAcceptance;
	static double acceptanceLimits[2] = {minimumAcceptance, maximumAcceptance};

	return acceptanceLimits;
}

