#include "TCanvas.h"
#include "TGraph.h"
#include "TSpline.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TH2D.h"
#include "TFile.h"
#include <vector>
#include "RegistrationEfficiency.h"
#include <TRandom3.h>
#include <TRandom.h>

using namespace std;
vector<vector<double>> getDetectionEfficiencyCorrectedEnergy(vector<vector<double>> );

int DetectionEfficiency() {

  vector<double> perEventWeight{};
  vector<vector<double>> perEventPhotonEnergies{};

  RegistrationEfficiency getPhasespaceEnergy;
  getPhasespaceEnergy.calculatePhasespaceEnergy(&perEventPhotonEnergies, &perEventWeight);

  vector<vector<double>> eventsAfterDetEfficiency = getDetectionEfficiencyCorrectedEnergy(perEventPhotonEnergies);

  return 0;

}


class DetectionEfficiencyInterpolator{
public:
	DetectionEfficiencyInterpolator(const char* prob_file = "regEffProbs.txt") {
		Graph = new TGraph( prob_file, "%lg %lg", ",");
		Spline_interpolator = new TSpline5("grs", Graph);
	}
	double getDetectionEfficiency(double inEnergy) { 
		return Spline_interpolator->Eval(inEnergy); //Provide energy in MeV
	}
private:
	TGraph *Graph = nullptr;
	TSpline5 *Spline_interpolator = nullptr;
};

// Calculates probabilities of gamma interaction in plastic scintillator based
// on gamma energy and saves it to a text file
void saveDetectionEfficiencyProbabilites(double detectorThickness = 2, const char *outFile = "regEffProbs.txt") {

	const int numberOfElements = 36;
	std::ofstream file(outFile);
	const double density = 1.03; // density of plastic
	double calculatedProbabilities[numberOfElements] = {0.0};

  double energies[numberOfElements] = {
      0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01,
      0.015, 0.02,   0.03,  0.04,  0.05,  0.06,  0.08,  0.1,   0.15,
      0.2,   0.3,    0.4,   0.5,   0.6,   0.8,   1,     1.25,  1.5,
      2,     3,      4,     5,     6,     8,     10,    15,    20}; //Energies in MeV

  double mass_attenuation_coeff[numberOfElements] = {
      2024,    640.9,   277,     82.7,   34.61,   17.53,   10.05,   4.22,
      2.204,   0.7705,  0.4358,  0.2647, 0.2194,  0.1997,  0.1881,  0.1736,
      0.1635,  0.1458,  0.1331,  0.1155, 0.1034,  0.09443, 0.08732, 0.07668,
      0.06894, 0.06166, 0.05611, 0.0481, 0.03848, 0.03282, 0.02907, 0.02641,
      0.0229,  0.02069, 0.0177,  0.01624};  //mass attenuation coefficients

  if (file.is_open())
  {
  for (int i = 0; i < numberOfElements; i++) {
	  calculatedProbabilities[i] = 1 - exp(-(mass_attenuation_coeff[i] * density * detectorThickness)); // probability of gamma interaction in scintillator
	  file << energies[i] << "," <<calculatedProbabilities[i] << "\n";
  } 
  file.close();
  }
}


vector<vector<double>> getDetectionEfficiencyCorrectedEnergy(vector<vector<double>> perEventPhotonEnergies)
{

  TRandom3* random = new TRandom3();
  random->SetSeed(0);

  const char *outFile = "regEffProbs.txt";
  const double detectorThickness = 2.0; // in cm
  saveDetectionEfficiencyProbabilites(detectorThickness, outFile);
  DetectionEfficiencyInterpolator detectionEfficiency(outFile);
 
  vector<vector<double>> eventsAfterDetEfficiency{};
  vector<vector<double>> perEventDetEfficiency{};

  const double detEffLowerLimit = 17.527;  
  const double detEffUpperLimit = 28.3146;

  TH2D* h1 = new TH2D("inEnergyVsDetEff", "Incoming Energy vs detection efficiency", 500, 0, 1, 100, 0, 100);
  auto pFile = new TFile("EnergyvsEFF.root ", "recreate");
 
  for(const auto& incomingEnergies : perEventPhotonEnergies)
  {
	  vector<double> energyAfterDetEfficiency{};
	  vector<double> perEnergyDetEfficiency{};
	  int count = 0;
	  double weight = 0.0;

	  for(const auto& incomingEnergy : incomingEnergies)
	  {
		  weight = 2;
		  if((incomingEnergy >= 0.105) && (incomingEnergy <= 0.511))
		  {
			  count++;
			  double calcDetEfficiency = detectionEfficiency.getDetectionEfficiency(incomingEnergy) * 100;
			 // perEnergyDetEfficiency.push_back(calcDetEfficiency);

			  while(true)
                          {
                                  double guessedDetEfficiency = random->Uniform(detEffLowerLimit,detEffUpperLimit);
                                  if(guessedDetEfficiency <= calcDetEfficiency)
                                  {
                                          h1->Fill(incomingEnergy, calcDetEfficiency);
                                          energyAfterDetEfficiency.push_back(incomingEnergy);
					  perEnergyDetEfficiency.push_back(calcDetEfficiency);
                                          break;
                                  }
			  }

		  }		  
	  
	  }
	  if(count == 4){
	  eventsAfterDetEfficiency.push_back(energyAfterDetEfficiency);
	  perEventDetEfficiency.push_back(perEnergyDetEfficiency);
	  perEnergyDetEfficiency.clear();
          energyAfterDetEfficiency.clear();
          count = 0;
	  }

	  else
	  {
		  perEnergyDetEfficiency.clear();
                  energyAfterDetEfficiency.clear();
		  count = 0;
	  }
	  
  }

  h1->GetXaxis()->SetTitle("Gamma_Icoming_Energy(MeV)");
  h1->GetYaxis()->SetTitle("Detection Efficiency #epsilon");
  h1->SetDirectory(gDirectory);
  pFile->Write();
  pFile->Close();
  delete pFile;

  return eventsAfterDetEfficiency;

}	
