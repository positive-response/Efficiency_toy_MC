#include "GeometricalEfficiency.cpp"
#include "RegistrationEfficiency.cpp"

int main()
{
double geometrical_eff = getGeometricalEfficieny();
double registration_eff = registrationEfficiency();

double total_eff = geometrical_eff * registration_eff;
std::cout<<"Total Efficiency:" << total_eff<<std::endl;

return 1;
}
