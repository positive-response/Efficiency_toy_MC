#include <iostream>
#include <vector>
#include <fstream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"

using namespace std;

double inter()
{
         double Energy[36] = {0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.25, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20};

        double mu[36] = {2024, 640.9, 277, 82.7, 34.61, 17.53, 10.05, 4.22, 2.204, 0.7705, 0.4358, 0.2647, 0.2194, 0.1997, 0.1881, 0.1736, 0.1635, 0.1458, 0.1331, 0.1155, 0.1034, 0.09443, 0.08732, 0.07668, 0.06894, 0.06166, 0.05611, 0.0481, 0.03848, 0.03282, 0.02907, 0.02641, 0.0229, 0.02069, 0.0177, 0.01624};

        double d = 20; //in cm
        double fun[36] = {0.0};
	double interpolated_mu = 0.0;

        for(int i = 0; i < 36; i++)
        {
         fun[i] = 1 - exp(-(mu[i]*d)); //probability of gamma interaction in plastic scintillator based on gamma energy
        }

//	TCanvas * c = new TCanvas("canvas", "E vs mu", 900, 800);
        TGraph* g = new TGraph(36, Energy, fun);
    //  g->Draw();
//	g->SetTitle("Gamma interaction probability; Energy(MeV); Probablity");
        TSpline5 *s = new TSpline5("grs",g); 
  //    s->SetLineColor(kRed);
//      s->Draw("same");
        interpolated_mu = s->Eval(0.511);
	return interpolated_mu;
}

