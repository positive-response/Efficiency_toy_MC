#include "TMath.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <iostream>
#include "TH3D.h"
#include "PhaseSpace.h"
#include "TFile.h"

using namespace std;

double getGeometricalEfficiency(const int, const int, vector<vector<double>>*, vector<double>*);
double* calculateAcceptance(const double, const double);

double GeometricalEfficiency(const int nParticles, const int minimumNoGamma){
	const double half_lengthOfDetector = 23.0;  //in cm
	const double radiusOfDetector = 42.5; // in cm
	double * thetaLimits = calculateAcceptance(radiusOfDetector, half_lengthOfDetector);
//	const double mininumTheta = *(thetaLimits + 0);
//	const double maximumTheta = *(thetaLimits + 1);

	std::vector<vector<double>> Energies;
	std::vector<double> EventWeight;
    
//	return getGeometricalEfficiency(mininumTheta, maximumTheta, nParticles, minimumNoGamma, &Energies);
return getGeometricalEfficiency(nParticles, minimumNoGamma, &Energies, &EventWeight);
}

/************************************************************************************************************/
double getGeometricalEfficiency( const int nDaughter = 5, const int minimumNumberOfGamma = 5, vector<vector<double>> *Energies = nullptr, vector<double> * EventWeight = nullptr)
{
	std::vector<vector<double>> gamma_theta_per_event;
	std::vector<vector<double>> gamma_phi_per_event;
	std::vector<double> perEventWeight;
	std::vector<vector<double>> perEventPhotonEnergies;	
    PhaseSpace getPhaseSpace4;
	getPhaseSpace4.numberOfDaughterParticles = nDaughter;
	getPhaseSpace4.getAnglesOfDaughterParticles(&gamma_theta_per_event, &gamma_phi_per_event, &perEventWeight);
	getPhaseSpace4.calculatePhasespaceEnergy(&perEventPhotonEnergies);  //energies are in MeV
        
	int total = 0;
	int accepted = 0;
	int n_min = 0;
	const double phi_min = -3.14;
	const double phi_max = 3.14;
	const double theta_min = 1.0472;
	const double theta_max = 2.0944;
	int r = 0;

	double gamma_theta = 0.0;
	double gamma_phi = 0.0;
	int rejected = 0;
	double wt = 0.0;

	vector<double> en{};
	for(int i = 0; i < static_cast<int>(perEventWeight.size()); i++)
	{
		gamma_theta = 0.0;
		gamma_phi = 0.0;
		n_min = 0;
        total++;
		r = 0;

		EventWeight->push_back(perEventWeight[i]);
		
		for(int j = 0; j < static_cast<int>(gamma_theta_per_event[i].size()); j++)
		{
			gamma_theta = gamma_theta_per_event[i][j];
		    gamma_phi = gamma_phi_per_event[i][j];
			
			if ((gamma_theta > theta_min) && (gamma_theta < theta_max) && (gamma_phi > phi_min) && (gamma_phi < phi_max))
			{
				wt = perEventWeight[i];
		   		en.push_back(perEventPhotonEnergies[i][j]);
			    n_min++;
			}
			else r++;
		}

		gamma_theta = 0.0;
		gamma_phi = 0.0; 
		if (n_min == minimumNumberOfGamma)
		{	accepted++;

			EventWeight->push_back(wt);
			Energies->push_back(en);
			wt = 0.0;
			en.clear();
		}
		else {
			rejected++;
			en.clear();
			wt = 0.0;
		}
		n_min = 0;
	}

//auto geo = new TFile("geo.root", "RECREATE");
//h2->SetDirectory(gDirectory);
//geo->Write();
//geo->Close();
//h2->Draw();
//h2->GetXaxis()->SetTitle("Angles (#theta) in radians");
//h2->GetYaxis()->SetTitle("counts");
cout<<"accepted: " << accepted<<" total: "<< total<<" rejected: "<< rejected<<endl;
std::cout<<"Geometrical Efficiency: " << static_cast<double>(accepted)/total<< std::endl;
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
