#include "TCanvas.h"
#include "TGraph.h"
#include "TSpline.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class DetectionEfficiencyInterpolator{
public:
	DetectionEfficiencyInterpolator(const char* prob_file = "regEffProbs.txt") {
		Graph = new TGraph( prob_file, "%lg %lg", ",");
		Spline_interpolator = new TSpline5("grs", Graph);
	}
	double getDetectionEfficiency(double inEnergy) { 
		return Spline_interpolator->Eval(inEnergy);
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
      2,     3,      4,     5,     6,     8,     10,    15,    20};

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

int main() {

  const char *outFile = "regEffProbs.txt";
  saveDetectionEfficiencyProbabilites(2.0, outFile); 
  DetectionEfficiencyInterpolator detectionEfficiency(outFile);
  std::cout << detectionEfficiency.getDetectionEfficiency(0.511) << std::endl;

  return 0;
}

