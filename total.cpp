#include "inter.cpp"
#include "Klein_4.cpp"
#include "phasespace4.cpp"

double total()
{
double geometrical_eff = phasespace4();
double registration_eff = Klein_4();
double detection_eff = inter();

double total_eff = geometrical_eff * registration_eff * detection_eff;

return total_eff;
}
