# Efficiency_toy_MC
Toy monte carlo to estimate total efficiency of J-PET Detector for positronium to 4 or 5 gamma decay.
Total efficiency is calculated as a product of registration efficiency, detection effieciency and geometrical efficiency.
Detection Efficiency:
The probability of gamma interacting with the plastic scintillator depends on the linear attenuation coefficient and the thickness of the plastic scintillator. The linear attenuation coefficient further depends on the energy of the incoming gamma. 
DetectionEfficiency.cpp calculates the detection efficiency of gamma based on this probability.
Geometrical Efficiency:
GeometricalEfficiency.cpp calculates the efficiency of jpet based on the acceptance of the detector.
Registration Efficiency:
RegistrationEfficiency.cpp calculates the registration efficiency of gamma from positronium decay taking into account the phasespace and the cross-section of gamma.
