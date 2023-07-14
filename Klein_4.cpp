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
using namespace std;

double klein(double , double );      // calculation of cross section
void phasespace4(vector<vector<double>> *, vector<double> *); // phasespace generation
double Efficiency(vector<double>*, const char* ); // efficiency calculation
vector<vector<double>> get_deposited_Energy(vector<vector<double>>& );
vector<double> get_small_deposited_energy(vector<vector<double>>& );

/********************Calculation of efficiency****************************/

double Efficiency(vector<double>* E_dep, const char* filename)
{
  const int nthr = 16;
  double Thr[nthr] = {0.0}; 
  int count[nthr] = {0};

  for(int k = 0; k < nthr; k++)
  {
	  Thr[k] = (k+1) * 360.0/nthr;
  }

  TEfficiency* pEff1 = 0;
  auto pFile = new TFile(filename, "recreate");
  pEff1 = new TEfficiency("eff1", "Threshold vs Efficiency; Threshold [in keV]; Efficiency #epsilon", 20, 0, 360);
  
  int c = 0;
  bool bPassed = false;
  double threshold = 0.0;

  for(int j = 0; j < nthr; j++)
  { 
	  c = 0;
	  bPassed = false;
	  threshold = Thr[j];
	 
	  for(double E : *E_dep)
	  {
		  bPassed = E > threshold;
		  if(bPassed){c++;}
		  pEff1->Fill(bPassed, threshold);
	  }

	  count[j] = c;
  }

  pEff1->SetDirectory(gDirectory);
  pFile->Write();
  pFile->Close();
  delete pFile;
/************************Efficiency corresponsing to 30 keV**************************/
  const double fixed_thr = 30.0; //in keV
  int ct = 0;
  double ratio = 0.0;
  int total_count1 = E_dep->size();
  for(auto En : *E_dep)
  {
	  if (fixed_thr > En)
		  ct++;
  }

  ratio = static_cast<double>(ct)/total_count1;
  std::cout<<"ratio:"<<ratio<<std::endl;
  std::cout<<"ct:"<<ct<<std::endl;
  std::cout<<"total_count1:"<<total_count1<<std::endl;
  return ratio; 
}

/**************************calculation of deposited energy********************************/
vector<vector<double>> get_deposited_Energy(vector<vector<double>>& events) { 
        TRandom3 * random = new TRandom3();
	random->SetSeed(0);

	double yguess = 0.0;
  	double theta = 0.0;
  	double yfun = 0.0;
	vector<vector<double>> events_with_deposited_e;
        vector<double> deposited_energies= {0,0,0,0};
	double scattered_e = 0;

	const double theta_max = 3.14159;
	const double theta_min = 0.1;
	const double cross_section_max = 0.489*1e-30;
	const double cross_section_min = 0.0;

  for(auto event : events)
  {
	  for(auto incoming_energy : event)
	  {
		  while(1){
                         theta = random->Uniform(theta_min, theta_max);
 			 yguess = random->Uniform(cross_section_min, cross_section_max);
			 yfun = klein(theta, incoming_energy);

			  if(yguess <= yfun)
			  {
				  scattered_e = incoming_energy/(1 + (1 - cos(TMath::RadToDeg()*theta)));
                    	          deposited_energies.push_back(incoming_energy -scattered_e);
                    		  break;
			  }
		
		  }
	  }
	  events_with_deposited_e.push_back(deposited_energies);
          deposited_energies.clear();
	  yguess = 0.0;
          theta = 0.0;
          yfun = 0.0;
  }
  return events_with_deposited_e;

}
/***********************calculation of smallest deposited energy**********************/
vector<double> get_small_deposited_energy(vector<vector<double>>& dep_energies)
{
	vector<double> small_deposited_energy;
	double small_edep = 0.0;

	for( auto Energies : dep_energies)
	{
		for(auto E_e : Energies)
		{
			small_edep = Energies.at(0);
			if (E_e  < small_edep)
                                small_edep = E_e;
		}
		small_deposited_energy.push_back(small_edep);
	}
	return small_deposited_energy;

}


/******************************************************************************/
double Klein_4()
{
	TRandom3 * random = new TRandom3();
	random->SetSeed(0);

       	TH1D* h3 = new TH1D("E_dep", "Energy", 500, 0, 550);
  	TH1D* h4 = new TH1D("E_dep_s", "EnergyS", 700, 0, 500);

  	vector<vector<double>> small_E_v;
  	vector<double> wts_v;
  	
  	const int nparticles = 4;
  	phasespace4(&small_E_v, &wts_v);
  	vector<vector<double>> Deposited_E = get_deposited_Energy(small_E_v);
	vector<double> Edep_small = get_small_deposited_energy(Deposited_E);
	vector<double> Smeared_E;

	const double smear_const = 0.044;
	double E = 0.0;
	double sigma = 0.0;
	double E_smear = 0.0;
	for(int i = 0; i < static_cast<int>(Edep_small.size()); i++)
	{
		E = Edep_small[i];
		sigma = (E * smear_const)/(sqrt(E/1000));
		E_smear = E + random->Gaus(0, sigma);
		Smeared_E.push_back(E_smear);
		h3->Fill(E, wts_v[i]);
		h4->Fill(E_smear, wts_v[i]);
	}

   double registration_eff = Efficiency(&Smeared_E, "smeared.root");
//   Efficiency(&Edep_small, "E_dep.root");
   h4->Draw();
   h3->SetTitle("Smallest energy deposition by gamma in 4gamma decay");
   h3->GetXaxis()->SetTitle("Energy_dep(keV)");
   h3->GetYaxis()->SetTitle("Counts");
   return registration_eff;
   }
/********************cross-section calculation**********************************/

double klein(double theta, double energy)
{
  const double pi =  3.14159;	
  const double N = 3.97 * 1e-30;    //(electron_radius^2)/2 (constant in the formula) in meters
  double alpha = energy / 511;     //energy of photon/rest mass energy of electron
  double cos_t = cos(theta);       
  double cos_one = 1 - cos_t;
  double cos_2 = 1 + pow(cos_t, 2);
  double KN = N * (cos_2 / (1 + (alpha*cos_one))) * (1 + (pow(cos_one, 2) * pow(alpha, 2)/ (cos_2 * (1 +(alpha * cos_one)))));
  double D_omega_by_d_theta = 2 * pi * sin(theta);
  double d_theta_by_d_E = pow((1 + alpha*(1 - cos_t)), 2)/(sin(theta)*alpha*energy);
  double result = d_theta_by_d_E * D_omega_by_d_theta * KN;
  return result;
}
/*****************************Phasespace calculation******************************************/


void phasespace4(vector<vector<double>> * small_E_v, vector<double> * wts_v)
{
 const double melectron = 0.000511; // mass of electron or positron in GeV
 TLorentzVector e(0.0,0.0,0.0,melectron); // four momenta of e-
 TLorentzVector p(0.0,0.0,0.0,melectron); //  four momenta of e+
 TLorentzVector W = e + p; //positronium atom

 const int npar = 4;
 double masses[npar] = {0.0};

 TGenPhaseSpace event;
 event.SetDecay(W, npar, masses);

 const int iter = 1000000;
 double weight = 0.0;
vector<double> gamma_E{};
 
 for(int i = 0; i < iter; i++)
   {
     weight = event.Generate();
     wts_v->push_back(weight);

     for(int j = 0; j< npar; j++)
       {
         TLorentzVector *gamma = event.GetDecay(j);
	 gamma_E.push_back(gamma->E() * 1000 * 1000); //From MeV to keV
       }
small_E_v->push_back(gamma_E); 
gamma_E.clear();
   }
}                                                                  
