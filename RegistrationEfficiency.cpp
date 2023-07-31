#include <iostream>
#include <TMath.h>
#include <cmath>
#include <algorithm>
#include <TGraph.h>
#include "TMath.h"
#include <TRandom.h>
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TFile.h"
#include "TEfficiency.h"
#include <TRandom3.h>
#include <vector>

#include "RegistrationEfficiency.h"

using namespace std;

double calculateCross_Section(double , double );      // calculation of cross section
//RegistrationEfficiency::calculatePhasespaceEnergy(vector<vector<double>> *, vector<double> *); // phasespace generation
double getRegistrationEfficiency(vector<double>*, const char*, const double ); // efficiency calculation
vector<vector<double>> getDepositedEnergy(vector<vector<double>>& );
vector<double> getSmallestDepositedEnergy(vector<vector<double>>& );


/********************Calculation of efficiency****************************/

double getRegistrationEfficiency(vector<double>* energiesWithSmearing, const char* filename, const double lowerThreshold = 30.0)
{
  const int numberOfThreshold = 20;
  double Thresholds[numberOfThreshold] = {0.0}; 

  for(int k = 0; k < numberOfThreshold; k++)
  {
	  Thresholds[k] = (k+1) * 360.0/numberOfThreshold;
  }

  TEfficiency* pEff1 = 0;
  auto pFile = new TFile(filename, "recreate");
  pEff1 = new TEfficiency("eff1", "Threshold vs Efficiency; Threshold [in keV]; Efficiency #epsilon", 20, 0, 360);
  
  bool bPassed = false;
  double threshold = 0.0;

  for(int j = 0; j < numberOfThreshold; j++)
  { 
	  bPassed = false;
	  threshold = Thresholds[j];
	 
	  for(double smearEnergy : *energiesWithSmearing)
	  {
		  bPassed = smearEnergy > threshold;
		  pEff1->Fill(bPassed, threshold);
	  }
  }

  pEff1->SetDirectory(gDirectory);
  pFile->Write();
  pFile->Close();
  delete pFile;
/************************Efficiency corresponsing to 30 keV**************************/
  int count = 0;
  double ratio = 0.0;
  int sizeOfEnergyVector = energiesWithSmearing->size();
  for(auto inEnergy : *energiesWithSmearing)
  {
	  if (inEnergy > lowerThreshold)
		  count++;
  }

  ratio = static_cast<double>(count)/sizeOfEnergyVector;
  count = 0;
  return ratio; 
}

/**************************calculation of deposited energy********************************/
vector<vector<double>> getDepositedEnergy(vector<vector<double>>& perEventPhotonEnergies) { 
        TRandom3 * random = new TRandom3();
	random->SetSeed(0);

	double guessedCross_section = 0.0;
  	double theta = 0.0;
  	double calculatedCross_section = 0.0;
	vector<vector<double>> perEventDepEnergy;
        vector<double> deposited_energies= {0,0,0,0};
	double scatteredGammaEnergy = 0;

	const double theta_max = 3.14159;
	const double theta_min = 0.1;
	const double maximumCross_Section = 0.489*1e-30;
	const double minimumCross_Section = 0.0;
	double alpha = 0.0;

  for(auto event : perEventPhotonEnergies)
  {
	  for(auto incoming_energy : event)
	  {
		  while(1){
			  alpha =  incoming_energy/0.511;
                          theta = random->Uniform(theta_min, theta_max);
 			  guessedCross_section = random->Uniform(minimumCross_Section, maximumCross_Section);
			  calculatedCross_section = calculateCross_Section(theta, incoming_energy);

	       		  if(guessedCross_section <= calculatedCross_section)
			  {
				  scatteredGammaEnergy = incoming_energy/(1 + alpha * (1 - cos(TMath::RadToDeg()*theta)));
				  deposited_energies.push_back(incoming_energy - scatteredGammaEnergy);
				  break;
			  }
		  }
	  }

	  perEventDepEnergy.push_back(deposited_energies);
          deposited_energies.clear();
	  guessedCross_section = 0.0;
          theta = 0.0;
          calculatedCross_section = 0.0;
  }
  return perEventDepEnergy;

}
/***********************calculation of smallest deposited energy**********************/
vector<double> getSmallestDepositedEnergy(vector<vector<double>>& perEvntDepositedEnergies)
{
	vector<double> smallestDepositedEnergy;
	double smallestEnergy = 0.0;

	for( auto perEventEnergies : perEvntDepositedEnergies)
	{
		for(auto  depositedGammaEnergy: perEventEnergies)
		{
			smallestEnergy = perEventEnergies.at(0);
			if (depositedGammaEnergy < smallestEnergy)
                                smallestEnergy = depositedGammaEnergy;
		}
		smallestDepositedEnergy.push_back(smallestEnergy * 1000);   //from MeV to keV
	}
	return smallestDepositedEnergy;

}

/******************************************************************************/
int main()
{
	TRandom3 * random = new TRandom3();
	random->SetSeed(0);

//    	TH1D* h3 = new TH1D("energiesWithSmearing", "Energy", 500, 0, 550);
//  	TH1D* h4 = new TH1D("energiesWithSmearing_s", "EnergyS", 700, 0, 500);

  	vector<vector<double>> perEventPhotonEnergies;
  	vector<double> perEventWeight;
	vector<double> energiesWithSmearing;
  	
  	RegistrationEfficiency::calculatePhasespaceEnergy(&perEventPhotonEnergies, &perEventWeight);  //energies are in MeV

  	vector<vector<double>> perEvntDepositedEnergies = getDepositedEnergy(perEventPhotonEnergies);
	vector<double> perEventSmallestDepositedEnergy = getSmallestDepositedEnergy(perEvntDepositedEnergies);
	

	const double smear_const = 0.044;
	for(int i = 0; i < static_cast<int>(perEventSmallestDepositedEnergy.size()); i++)
	{
		double energy = perEventSmallestDepositedEnergy[i];
		double sigma = (energy * smear_const)/(sqrt(energy/1000));
		double smearedEnergy = energy + random->Gaus(0, sigma);

		energiesWithSmearing.push_back(smearedEnergy);
	//	h3->Fill(E, perEventWeight[i]);
	//	h4->Fill(E_smear, perEventWeight[i]);
	}
	double registration_eff = getRegistrationEfficiency(&energiesWithSmearing, "smeared.root", 30.0);
   //   h4->Draw();
   //   h3->SetTitle("Smallest energy deposition by gamma in 4gamma decay");
   //   h3->GetXaxis()->SetTitle("Energy_dep(keV)");
   //   h3->GetYaxis()->SetTitle("Counts");
        std::cout<<"registration_efficiency: "<< registration_eff<<std::endl;
        return 0;

   }
/********************cross-section calculation**********************************/

double calculateCross_Section(double theta, double energy)
{
  const double pi =  3.14159;	
  const double numericalConst = 3.97 * 1e-30;    //(electron_radius^2)/2 (constant in the formula) in meters
  double alpha = energy / 511;     //energy of photon/rest mass energy of electron
  double cosTheta = cos(theta);       
  double oneMinusCosTheta = 1 - cosTheta;
  double cos_2 = 1 + pow(cosTheta, 2);
  double dSigma_by_dOmega = numericalConst * (cos_2 / (1 + (alpha*oneMinusCosTheta))) * (1 + (pow(oneMinusCosTheta, 2) * pow(alpha, 2)/ (cos_2 * (1 +(alpha * oneMinusCosTheta)))));
  double dOmega_by_dTheta = 2 * pi * sin(theta);
  double dTheta_by_dE = pow((1 + alpha*(1 - cosTheta)), 2)/(sin(theta)*alpha*energy);
  double differentialCross_section =  dSigma_by_dOmega * dOmega_by_dTheta * dTheta_by_dE;
  return differentialCross_section;
}
/*****************************Phasespace calculation******************************************/
/*

int RegistrationEfficiency::calculatePhasespaceEnergy(std::vector<std::vector<double>> * perEventPhotonEnergies, std::vector<double> * perEventWeight)
{
 const double electron_mass = 0.000511; // mass of electron or positron in GeV
 TLorentzVector electron(0.0,0.0,0.0,electron_mass); // four momenta of e-
 TLorentzVector proton(0.0,0.0,0.0,electron_mass); //  four momenta of e+
 TLorentzVector parent_particle = electron + proton; //positronium atom

 const int numberOfDaughterParticles = 4;
 double masses[numberOfDaughterParticles] = {0.0};

 TGenPhaseSpace event;
 event.SetDecay(parent_particle, numberOfDaughterParticles, masses);

 const int iteration = 1000000;
 double weight = 0.0;
 vector<double>energiesOfGamma{};
 
 for(int i = 0; i < iteration; i++)
   {
     weight = event.Generate();
     perEventWeight->push_back(weight);

     for(int j = 0; j < numberOfDaughterParticles; j++)
       {
         TLorentzVector *gamma = event.GetDecay(j);
	 energiesOfGamma.push_back(gamma->E() * 1000 * 1000); //From MeV to keV
       }
perEventPhotonEnergies->push_back(energiesOfGamma); 
energiesOfGamma.clear();
   }
   return 0;
}*/
