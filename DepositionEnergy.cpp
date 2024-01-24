#include <TMath.h>
#include <cmath>
#include <algorithm>
#include <TGraph.h>
#include "TMath.h"
#include <TRandom.h>
#include "TH1D.h"
#include "TFile.h"
#include <TRandom3.h>
#include <vector>

using namespace std;

/**********************************************************************/
double DepositionEnergy()
{
	const int nthr = 45;
	const double upperEnergy = 550.0;
	double energies[nthr] = {0.0};
	double depositedEnergies[nthr] = {0.0};
	for(int i = 0; i < nthr; i++)
	{
		energies[i] = (i+1) * upperEnergy/nthr;
	}

	double alpha = 0.0;
	const double theta = 180;

	for(int j = 0; j < nthr; j++)
	{
		alpha = energies[j]/511;
		depositedEnergies[j] = energies[j] - (energies[j]/(1 +  (alpha* (1 - std::cos(theta*TMath::DegToRad())))));

	//	std::cout<<"dep: "<< depositedEnergies[j]<<" for energy: "<<energies[j]<<std::endl;
		alpha = 0.0;
	}
	TGraph* g = new TGraph(nthr, energies, depositedEnergies);
	g->GetXaxis()->SetTitle("Incoming gamma Energy(keV)");
	g->GetYaxis()->SetTitle("Deposited Energy(keV)");
	g->SetTitle("Incoming gamma energy vs Energy deposited by gamma");

	g->Draw("ac");
	return cos(theta*TMath::DegToRad());
	}
