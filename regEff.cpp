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
#include "DetectionEfficiency.cpp"

using namespace std;

double calculateCross_Section(double , double );      // calculation of cross section
//RegistrationEfficiency::calculatePhasespaceEnergy(vector<vector<double>> *, vector<double> *); // phasespace generation
double getRegistrationEfficiency(vector<double>*, const char*, const double ); // efficiency calculation
vector<vector<double>> getDepositedEnergy(vector<vector<double>>& );
vector<double> getSmallestDepositedEnergy(vector<vector<double>>& );

TH1D* functionToGetSmallestEnergy(const int no)
{
	TRandom3 * random = new TRandom3();
        random->SetSeed(0);
	TH1D* h3 = new TH1D("smallestEnergiesWithSmearingfor4gamma", "EnergyS", 200, 0, 200);

	vector<vector<double>> perEventPhotonEnergies;
        vector<double> perEventWeight;
        vector<double> energiesWithSmearing;
	const int npar = no;
	RegistrationEfficiency getPhasespaceEnergies(0.000511, npar ,1000000);

	getPhasespaceEnergies.calculatePhasespaceEnergy(&perEventPhotonEnergies, &perEventWeight);  //energies are in MeV
	vector<vector<double>> detEfficiencyCorrectedEvents = getDetectionEfficiencyCorrectedEnergy(perEventPhotonEnergies, npar);
  	vector<vector<double>> perEvntDepositedEnergies = getDepositedEnergy(detEfficiencyCorrectedEvents);
	vector<double> perEventSmallestDepositedEnergy = getSmallestDepositedEnergy(perEvntDepositedEnergies);

        const double smear_const = 0.044;
        const double registrationThreshold = 30.0;
        for(int i = 0; i < static_cast<int>(perEventSmallestDepositedEnergy.size()); i++)
        {
                double energy = perEventSmallestDepositedEnergy[i];
                double sigma = (energy * smear_const)/(sqrt(energy/1000));
                double smearedEnergy = energy + random->Gaus(0, sigma);

                energiesWithSmearing.push_back(smearedEnergy);
                h3->Fill(smearedEnergy, perEventWeight[i]);
        }
        double registration_eff = getRegistrationEfficiency(&energiesWithSmearing, "smeared.root", registrationThreshold);
//      double detectionEfficiency = static_cast<double>(detEfficiencyCorrectedEvents.size())/perEventPhotonEnergies.size();
       h3->SetTitle("Smallest energy deposition by gamma in o-Ps decay(Detection Efficiency included)");
	h3->GetXaxis()->SetTitle("Deposited Energy [keV]");
        h3->GetYaxis()->SetTitle("Counts");
	return h3;
}

double regEff()
{

        TH1D* gamma4 = functionToGetSmallestEnergy(4);
	TH1D* gamma5 = functionToGetSmallestEnergy(5);

        TCanvas* c2 = new TCanvas("c2", "new canvas");

        TH1D* h1C = (TH1D*)gamma4->Clone("h1C");
        gamma4->Scale(1000 / h1C->Integral());
        gamma4->Draw();

	gamma5->Draw("SAME");
	gamma4->SetLineColor(kRed);
	gamma4->SetLineWidth(2);

	gamma5->SetLineColor(kBlue);
	gamma5->SetLineWidth(2);


	auto legend = new TLegend(0.1, 0.3, 0.3, 0.5);
	legend->AddEntry(gamma4, "Smallest deposited energy for 4 #gamma");
	legend->AddEntry(gamma5, "Smallest deposited energy for 5 #gamma");
	legend->Draw();
        return 0.0;

//     h3->SetTitle("Smallest energy deposition by gamma in 4gamma decay(Detection Efficiency included)");
//      h4->GetXaxis()->SetTitle("Deposited Energy [keV]");
//      h4->GetYaxis()->SetTitle("Counts");
//        std::cout<<"Registration_efficiency(after Detection Efficiency corrected): "<< registration_eff<<std::endl;
//      std::cout<<"Fraction of events within the allowed probability range(Detection Efficiency): " << detectionEfficiency<<std::endl;	
//      return registration_eff;

   }
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
				  scatteredGammaEnergy = incoming_energy/(1 + alpha * (1 - cos(theta)));
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

