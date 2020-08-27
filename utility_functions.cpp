// implementation of utility functions class

#include "utility_functions.h"

#include <cmath>

// destructor
utility_functions::~utility_functions(){
	//if fScintillation_function != nullptr {
		delete fScintillation_function_muon;
		delete fScintillation_function_electron;
		delete fScintillation_function_alpha;
	//}
}

// poisson distribution function
int utility_functions::poisson(const double mean, const double draw, const double eng) const {
	int number = 0;
	const int border = 16;
	double limit = 2e9;

	if(mean <= border) {
		double position = draw;
		double poissonValue = std::exp(-mean);
		double poissonSum = poissonValue;

		while(poissonSum <= position) {
			number++;
			poissonValue *= mean/number;
			poissonSum += poissonValue;
		}
		return number;
	} // the case of mean <= 16

	double value, t, y;
	t = std::sqrt(-2*std::log(draw));
	y = 2*3.141592654*eng;
	t *= std::cos(y);
	value = mean + t*std::sqrt(mean) + 0.5;
	if(value <= 0) {return 0; }
	if(value >= limit) { return limit; }
	return value;
}

// scintillation time spectrum function
double utility_functions::scintillation_function(const double *t, const double *par) {
	double time = *t;
	double t_singlet = par[0];
	double t_triplet = par[1];
	double type = par[2]; 		// type will be defined as 0,1 or 2:  0 is a muon, 1 is an electron, 2 is an alpha particle
	double singlet_part;
	double triplet_part;

	if(type == 0){ // particle is a muon
	  singlet_part = 0.23;
	  triplet_part = 0.77;
	}
	if(type == 1){ // particle is an electron
	  singlet_part = 0.27;
	  triplet_part = 0.73;
	}
	if(type == 2){ // particle is an alpha
	  singlet_part = 0.70;
	  triplet_part = 0.30;
	}

	double Scintillation = exp(-(time/t_singlet))*singlet_part/t_singlet + exp(-(time/t_triplet))*triplet_part/t_triplet;

	return Scintillation;
}

// function to create scintillation function TF1s with required parameters
void utility_functions::initalise_scintillation_functions(const double t_singlet, const double t_triplet, const double scint_time_window) {

	// create scintillation spectrum
	// muon
	fScintillation_function_muon = new TF1("Scintillation Timing Muon", scintillation_function, 0, scint_time_window, 3);
	fScintillation_function_muon->SetParameter(0, t_singlet); 
    fScintillation_function_muon->SetParameter(1, t_triplet);  
    fScintillation_function_muon->FixParameter(2, 0); 	// muon
    fScintillation_function_muon->SetNpx(50000);     
    // electron
    fScintillation_function_electron = new TF1("Scintillation Timing Electron", scintillation_function, 0, scint_time_window, 3);
	fScintillation_function_electron->SetParameter(0, t_singlet); 
    fScintillation_function_electron->SetParameter(1, t_triplet);  
    fScintillation_function_electron->FixParameter(2, 1); 	// electron
    fScintillation_function_electron->SetNpx(50000);     
    // alpha
    fScintillation_function_alpha = new TF1("Scintillation Timing Alpha", scintillation_function, 0, scint_time_window, 3);
	fScintillation_function_alpha->SetParameter(0, t_singlet); 
    fScintillation_function_alpha->SetParameter(1, t_triplet);  
    fScintillation_function_alpha->FixParameter(2, 2); 	// alpha 
    fScintillation_function_alpha->SetNpx(50000);    
}

double utility_functions::get_scintillation_time(const int &particle_type) {
	// sample function based on particle type
	if (particle_type == 0) return fScintillation_function_muon->GetRandom();
	else if (particle_type == 1) return fScintillation_function_electron->GetRandom();
	else if (particle_type == 2) return fScintillation_function_alpha->GetRandom();
	else return 0;
}