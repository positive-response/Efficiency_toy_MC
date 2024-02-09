#ifndef PHASESPACE_H
#define PHASESPACE_H

#include <vector>
#include <TMath.h>
#include <cmath>
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

using namespace std;
class PhaseSpace
{
    public:
	    int iteration = 1000000;
	    double electron_mass = 0.000511;
	    double numberOfDaughterParticles = 5;
	    
	    PhaseSpace()
		    {
			    numberOfDaughterParticles = 5;
			    iteration = 1000000;
			    electron_mass = 0.000511;

		    }
		    
	    TGenPhaseSpace generateEvents(const int numberOfDaughterParticles)
	    {
		    TGenPhaseSpace event;
    	            TLorentzVector electron(0.0,0.0,0.0,electron_mass); // four momenta of e-
		    TLorentzVector proton(0.0,0.0,0.0,electron_mass); //  four momenta of e+
		    TLorentzVector parent_particle = electron + proton; //positronium atom
		    std::vector<double> masses(numberOfDaughterParticles, 0.0);
		    event.SetDecay(parent_particle, numberOfDaughterParticles, masses.data());
		    return event;
	    }
	    void calculatePhasespaceEnergy(std::vector<std::vector<double>> * perEventPhotonEnergies/*, std::vector<double> * perEventWeight*/)
	    {
		    TGenPhaseSpace event = generateEvents(numberOfDaughterParticles);
		    double weight = 0.0;
 	            vector<double>energiesOfGamma{};

		for(int i = 0; i < iteration; i++)
		{
			weight = event.Generate();
//			perEventWeight->push_back(weight);
			for(int j = 0; j < numberOfDaughterParticles ; j++)
			{
				TLorentzVector *gamma = event.GetDecay(j);
				energiesOfGamma.push_back(gamma->E() * 1000); //From GeV to MeV
			}
			perEventPhotonEnergies->push_back(energiesOfGamma); 
			energiesOfGamma.clear();
		}
	    }

    
	    void getAnglesOfDaughterParticles(std::vector<std::vector<double>> * gamma_theta_per_event,  std::vector<std::vector<double>> *gamma_phi_per_event ,std::vector<double> * perEventWeight)
	    {
		    TGenPhaseSpace event = generateEvents(numberOfDaughterParticles);
		    double weight = 0.0;
		    vector<double>  gamma_theta{};
		    vector<double>  gamma_phi{};
		    for(int i = 0; i < iteration; i++)
		    {
			    weight = event.Generate();
			    perEventWeight->push_back(weight);
		    
			    for(int j = 0; j< numberOfDaughterParticles; j++)
			    {	    
				    TLorentzVector *gamma = event.GetDecay(j);
				    gamma_theta.push_back(gamma->Theta());
				    gamma_phi.push_back(gamma->Phi());
			    }
			    gamma_theta_per_event->push_back(gamma_theta);
			    gamma_phi_per_event->push_back(gamma_phi);
			    gamma_theta.clear();
			    gamma_phi.clear();
		    }
	    }

		
};

#endif
