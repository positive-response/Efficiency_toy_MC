# Efficiency_toy_MC
Toy monte carlo to estimate total efficiency of J-PET Detector for positronium to 4 gamma decay.
Detection Efficiency:
The probability of gamma interacting with the plastic scintillator depends on the linear attenuation coefficient and the thickness of the plastic scintillator. The linear attenuation coefficient further depends on the energy of the incoming gamma. 
inter.cpp calculates the detection efficiency of gamma based on this probability.
Geometrical Efficiency:
phasespace4.cpp calculates the efficiency of jpet based on the acceptance of the detector.
Registration Efficiency:
Klein_4.cpp calculates the registration efficiency of gamma from positronium decay taking into account the phasespace and the cross-section of gamma.
