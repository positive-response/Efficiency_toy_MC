#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TH2D.h"
#include "TFile.h"
#include <vector>
#include "PhaseSpace.h"
#include "DetectionEfficiencyInterpolator.h"

using namespace std;
vector<vector<double>> getDetectionEfficiencyCorrectedEnergy(vector<vector<double>> , const int);
/*
int DetectionEfficiency()
{
	PhaseSpace phase;
	phase.numberOfDaughterParticles = 5;
	phase.iteration = 20;
	vector<vector<double>> perEventPhotonEnergies;
	vector<double>perEventWeight;

	phase.calculatePhasespaceEnergy(&perEventPhotonEnergies, &perEventWeight);
	
	getDetectionEfficiencyCorrectedEnergy(perEventPhotonEnergies, 1);
	return 1;
}
*/
vector<vector<double>> getDetectionEfficiencyCorrectedEnergy(vector<vector<double>> perEventPhotonEnergies, const int nparticlesDetected)
{

  TRandom3* random = new TRandom3();
  random->SetSeed(0);

  const char *outFile = "regEffProbs.txt";
  const double detectorThickness = 2.0; // in cm
  saveDetectionEfficiencyProbabilites(detectorThickness, outFile);
  DetectionEfficiencyInterpolator detectionEfficiency(outFile);
 
  std::vector<std::vector<double>> eventsAfterDetEfficiency{};
  std::vector<std::vector<double>> perEventDetEfficiency{};

//  const double detEffLowerLimit = 17.527; //efficiency corresponding to 105 kev 
//  const double detEffUpperLimit = 28.3146; //efficiency corresponding to 511 kev

  const double detEffLowerLimit = 0.000165; //efficiency corresponding to 105 kev 
  const double detEffUpperLimit = 0.0018; //efficiency corresponding to 511 kev

  TH2D* h1 = new TH2D("inEnergyVsDetEff", "Incoming Energy vs detection efficiency", 500, 0, 1, 100, 0, 100);
  auto pFile = new TFile("EnergyvsEFF.root ", "recreate");
 
  int Ncount = 0;
  int count = 0;
  for(const auto& incomingEnergies : perEventPhotonEnergies)
  {
	  count++;
	  vector<double> energyAfterDetEfficiency{};
	  vector<double> perEnergyDetEfficiency{};
	   vector<double> efficiencyProd{};
	   double weight = 1.0;
	  for(const auto& incomingEnergy : incomingEnergies)
	  {	  
		  double calcDetEfficiency = detectionEfficiency.getDetectionEfficiency(incomingEnergy);
		  weight*=calcDetEfficiency;
		  h1->Fill(incomingEnergy, calcDetEfficiency);
		  energyAfterDetEfficiency.push_back(incomingEnergy);
		  perEnergyDetEfficiency.push_back(calcDetEfficiency);
	  }
			  double guessedDetEfficiency = random->Uniform(0, 1);
			  if(guessedDetEfficiency < weight)
			  {
				  Ncount++;
				  eventsAfterDetEfficiency.push_back(energyAfterDetEfficiency);
				  perEventDetEfficiency.push_back(perEnergyDetEfficiency);

				  perEnergyDetEfficiency.clear();
				  energyAfterDetEfficiency.clear();
				  weight = 1.0;
			  }	  
  
	  else{
		  perEnergyDetEfficiency.clear();
		  energyAfterDetEfficiency.clear();
		  weight = 1.0;
		  
	  }	 	  
  }	

  std::cout<<"Event count after det condition:" << Ncount <<endl;
  std::cout<<"Event count before det condition:" << count  <<endl;
  h1->GetXaxis()->SetTitle("Gamma_Icoming_Energy(MeV)");
  h1->GetYaxis()->SetTitle("Detection Efficiency #epsilon");
  h1->SetDirectory(gDirectory);
  pFile->Write();
  pFile->Close();
  delete pFile;
  return eventsAfterDetEfficiency;

}	
