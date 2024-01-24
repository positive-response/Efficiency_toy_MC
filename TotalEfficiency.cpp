#include "GeometricalEfficiency.cpp"
#include "RegistrationEfficiency.cpp"

int main(int argc, char* argv[])
{
	const int decayMode = atoi(argv[1]);
        const int nRequiredDecayParticles = atoi(argv[2]);
	std::cout<<"***********************stage 1 started: ***********************"<<std::endl;
	double geometrical_eff = GeometricalEfficiency(decayMode, nRequiredDecayParticles);
	std::cout<<"********************stage 2 started:*********************** "<<endl;
	double registration_eff = RegistrationEfficiency(decayMode, decayMode);
	double total_eff = geometrical_eff * registration_eff;
	std::cout<<"Total Efficiency(in %):" << total_eff*100<<std::endl;
	std::cout<<"**********************calculation finished********************************"<<std::endl;
	return 1;
}
