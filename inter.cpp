#include "TCanvas.h"
#include "TGraph.h"
#include "TSpline.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class RegistrationEfficiencyInterpolator {
public:
  RegistrationEfficiencyInterpolator(double *energies, double *probValues) {
    g = new TGraph(numberOfElements, energies, probValues);
    s = new TSpline5("grs", g);
  }

  double getRegistrationEfficiency(double inEnergy) { return -1; }

private:
  const int numberOfElements = 36;
  TGraph *g = nullptr;
  TSpline5 *s = nullptr;
};

// Calculates probabilities of gamma interaction in plastic scintillator based
// on gamma energy and saves it to a text file
void saveRegistrationEfficiencyProbabilites(
    double detectorThickness = 2, const char *outFile = "regEffProbs.txt") {
  const int numberOfElements = 36;
  double Energy[numberOfElements] = {
      0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01,
      0.015, 0.02,   0.03,  0.04,  0.05,  0.06,  0.08,  0.1,   0.15,
      0.2,   0.3,    0.4,   0.5,   0.6,   0.8,   1,     1.25,  1.5,
      2,     3,      4,     5,     6,     8,     10,    15,    20};

  double mu[numberOfElements] = {
      2024,    640.9,   277,     82.7,   34.61,   17.53,   10.05,   4.22,
      2.204,   0.7705,  0.4358,  0.2647, 0.2194,  0.1997,  0.1881,  0.1736,
      0.1635,  0.1458,  0.1331,  0.1155, 0.1034,  0.09443, 0.08732, 0.07668,
      0.06894, 0.06166, 0.05611, 0.0481, 0.03848, 0.03282, 0.02907, 0.02641,
      0.0229,  0.02069, 0.0177,  0.01624};

  const double density = 1.03; // density of plastic
  double fun[numberOfElements] = {0.0};
  double interpolated_mu = 0.0;

  for (int i = 0; i < numberOfElements; i++) {
    fun[i] =
        1 -
        exp(-(mu[i] * density *
              detectorThickness)); // probability of gamma interaction in
                                   // plastic scintillator based on gamma energy
  }
  /// Save Energies and fun to file
}

// double inter1(double en)
//{

////	TCanvas * c = new TCanvas("canvas", "E vs mu", 900, 800);
// TGraph* g = new TGraph(36, Energy, fun);
////  g->Draw();
////	g->SetTitle("Gamma interaction probability; Energy(MeV); Probablity");
// TSpline5 *s = new TSpline5("grs",g);
////    s->SetLineColor(kRed);
////      s->Draw("same");
// interpolated_mu = s->Eval(en);
// return interpolated_mu;
//}

int main() {
  const int numberOfElements = 36;
  // std::vector<double> energies = {};
  // energies.size() ==36;
  // RegistrationEfficiencyInterpolator (const std::vector<double>& energies,
  // ..)
  double energies[numberOfElements] = {
      0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.01,
      0.015, 0.02,   0.03,  0.04,  0.05,  0.06,  0.08,  0.1,   0.15,
      0.2,   0.3,    0.4,   0.5,   0.6,   0.8,   1,     1.25,  1.5,
      2,     3,      4,     5,     6,     8,     10,    15,    20};
  double FAKEPROBABILITIES[numberOfElements] = {
      2024,    640.9,   277,     82.7,   34.61,   17.53,   10.05,   4.22,
      2.204,   0.7705,  0.4358,  0.2647, 0.2194,  0.1997,  0.1881,  0.1736,
      0.1635,  0.1458,  0.1331,  0.1155, 0.1034,  0.09443, 0.08732, 0.07668,
      0.06894, 0.06166, 0.05611, 0.0481, 0.03848, 0.03282, 0.02907, 0.02641,
      0.0229,  0.02069, 0.0177,  0.01624};

  RegistrationEfficiencyInterpolator interpolator(energies, FAKEPROBABILITIES);
  std::cout << interpolator.getRegistrationEfficiency(0.511) << std::endl;
  // interpolator.getRegistrationEfficiency(0.511);

  //	double en_511 = inter1(0.511);
  //	double en_1274 = inter1(1.274);
  // std::cout<<"511: " <<inter1(0.511)<<"  "<<"1274:
  // "<<inter1(1.274)<<std::endl;
  return 0;
}
