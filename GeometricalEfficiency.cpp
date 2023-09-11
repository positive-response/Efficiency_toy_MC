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
        std::cout << getGeometricalEfficieny(mininumTheta, maximumTheta, minimumNumberOfGamma) << std::endl;
	return 0;
}

/************************************************************************************************************/
double getGeometricalEfficieny(const double theta_min = 1.0472, const double theta_max = 2.0944, const int minimumNumberOfGamma = 3)
{
TH3D* h1 = new TH3D("h1", "Position of gamma before geometrical cut", 100, -1, 1, 100,-1,1, 100, -1, 1);
TH3D* h2 = new TH3D("h2", "position of gamma after geometrical cut", 100, -1, 1, 100, -1, 1, 100,-1,1);
TH1D* h3 = new TH1D("h3", "Largest angle distribution", 200, 0, 3.50);

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
	   double small_theta = gamma_theta[0];
	   for(int j = 0; j< numberOfDaughterParticles; j++)
	   {
		   TLorentzVector *gamma = event.GetDecay(j);
                   gamma_theta[j] = gamma->Theta();
		   h1->Fill(cos(gamma->Theta())*sin(gamma->Phi()), sin(gamma->Theta())*sin(gamma->Phi()),cos(gamma->Phi()), weight);
                     //  h1->Fill(gamma->Z(), gamma->Y(), weight);
		
		  if ((gamma_theta[j] > theta_min) & (gamma_theta[j] < theta_max))
		  {	  n_min++;
	
		  h2->Fill(cos(gamma->Theta())*sin(gamma->Phi()), sin(gamma->Theta())*sin(gamma->Phi()),cos(gamma->Phi()), weight);
		  }
	   }

	   small_theta = gamma_theta[0];
	     for (int x = 1; x < numberOfDaughterParticles; x++)
		         {
				 if (gamma_theta[x] > small_theta)
				 {
					 small_theta = gamma_theta[x];
				 }
			 }
	     h3->Fill(small_theta, weight);

	     if (n_min >= minimumNumberOfGamma)
		     accepted++;
	     n_min = 0;
  }
h3->Draw();

h3->GetXaxis()->SetTitle("Angles (#theta) in radians");
h3->GetYaxis()->SetTitle("counts");
std::cout<<"Geometrical Efficiency: " << static_cast<double>(accepted)/total<< std::endl;
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

