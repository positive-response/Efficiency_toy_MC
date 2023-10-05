#include "TMath.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <iostream>
#include "TH3D.h"

using namespace std;

double getGeometricalEfficieny(const double , const double, const int);
double* calculateAcceptance(const double, const double);

int GeometricalEfficiency(){
	const double half_lengthOfDetector = 25.0;  //in cm
	const double radiusOfDetector = 42.5; // in cm
	double * thetaLimits = calculateAcceptance(radiusOfDetector, half_lengthOfDetector);
	const double mininumTheta = *(thetaLimits + 0);
	const double maximumTheta = *(thetaLimits + 1);
	const int minimumNumberOfGamma = 3;
       std::cout << getGeometricalEfficieny(mininumTheta, maximumTheta, minimumNumberOfGamma)<<std::endl;

      // std::cout << getGeometricalEfficieny(mininumTheta, maximumTheta, 2) << " 2gamma" <<std::endl;
      // std::cout << getGeometricalEfficieny(mininumTheta, maximumTheta, 3) <<" 3gamma" << std::endl;
      // std::cout << getGeometricalEfficieny(mininumTheta, maximumTheta, 4) <<" 4gamma" << std::endl;
	delete[] thetaLimits;
	return 0;

}

/************************************************************************************************************/
double getGeometricalEfficieny(const double theta_min = 1.0472, const double theta_max = 2.0944, const int minimumNumberOfGamma = 4)
{
//TH3D* h1 = new TH3D("h1", "Position of gamma before geometrical cut", 100, -1, 1, 100,-1,1, 100, -1, 1);
TH3D* h2 = new TH3D("h2", "position of gamma after geometrical cut", 100, -1, 1, 100, -1, 1, 100,-1,1);
//TH1D* h3 = new TH1D("h3", "Largest angle distribution", 200, 0, 3.50);

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
double gamma_phi[numberOfDaughterParticles]={0.0};
const double phi_min = -3.14;
const double phi_max = 3.14;
int total = 0;
int accepted = 0;
int n_min = 0;
for(int i = 0; i < iteration; i++)
   {
	   n_min = 0;
	   total++;
	   weight = event.Generate();
	   double small_theta = gamma_theta[0];
	   for(int j = 0; j< numberOfDaughterParticles; j++)
	   {
		   TLorentzVector *gamma = event.GetDecay(j);
                   gamma_theta[j] = gamma->Theta();
		   gamma_phi[j] = gamma->Phi();
//		   h1->Fill(cos(gamma->Theta())*sin(gamma->Phi()), sin(gamma->Theta())*sin(gamma->Phi()),cos(gamma->Phi()), weight);
		
		  if ((gamma_theta[j] > theta_min) & (gamma_theta[j] < theta_max) & (gamma_phi[j] > phi_min) & (gamma_phi[j] < phi_max))
		  {	  n_min++;
	
		  h2->Fill(cos(gamma->Theta())*sin(gamma->Phi()), sin(gamma->Theta())*sin(gamma->Phi()),cos(gamma->Phi()), weight);
		  }
	   }
/*
	   small_theta = gamma_theta[0];
	     for (int x = 1; x < numberOfDaughterParticles; x++)
		         {
				 if (gamma_theta[x] > small_theta)
				 {
					 small_theta = gamma_theta[x];
				 }
			 }
	     h3->Fill(small_theta, weight);
*/
	     if (n_min >= minimumNumberOfGamma)
		     accepted++;
	     n_min = 0;
  }
h2->Draw();
h2->GetXaxis()->SetTitle("Angles (#theta) in radians");
h2->GetYaxis()->SetTitle("counts");
std::cout<<"Geometrical Efficiency: " << static_cast<double>(accepted)/total<< std::endl;
//delete h2;
return static_cast<double>(accepted)/total;
}
/*****************************************************************************************/
double* calculateAcceptance(const double radius = 42.5, const double half_length = 25.0)
{

	const double minimumAcceptanceTheta = TMath::ATan(radius/half_length);
	const double maximumAcceptanceTheta = 2 * minimumAcceptanceTheta; 
	static double acceptanceLimits[2] = {minimumAcceptanceTheta, maximumAcceptanceTheta};

	return acceptanceLimits;
}

