#ifndef REGISTRATIONEFFICIENCY_H
#define REGISTRATIONEFFICIENCY_H

#include <vector>
#include <TMath.h>
#include <cmath>
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

using namespace std;
class RegistrationEfficiency
{
    public:
    RegistrationEfficiency(const double eMass = 0.000511, const double nParticles = 4, const int iter = 1000000){

    electron_mass = eMass;
    numberOfDaughterParticles = nParticles;
    iteration = iter;
    };

    int calculatePhasespaceEnergy(std::vector<std::vector<double>> * perEventPhotonEnergies, std::vector<double> * perEventWeight)
     {
       	TLorentzVector electron(0.0,0.0,0.0,electron_mass); // four momenta of e-
	TLorentzVector proton(0.0,0.0,0.0,electron_mass); //  four momenta of e+
	TLorentzVector parent_particle = electron + proton; //positronium atom
	double masses[numberOfDaughterParticles] = {0.0};
	TGenPhaseSpace event;
	event.SetDecay(parent_particle, numberOfDaughterParticles, masses);

        double weight = 0.0;
 	vector<double>energiesOfGamma{};
 
 	for(int i = 0; i < iteration; i++)
	{
		weight = event.Generate();
     		perEventWeight->push_back(weight);
		for(int j = 0; j < numberOfDaughterParticles ; j++)
		{
			TLorentzVector *gamma = event.GetDecay(j);
			energiesOfGamma.push_back(gamma->E() * 1000); //From GeV to MeV
		}
		perEventPhotonEnergies->push_back(energiesOfGamma); 
		energiesOfGamma.clear();
	}
	return 0;
     }
    private:

    double electron_mass;
    int numberOfDaughterParticles;
    int iteration;


};
#endif
