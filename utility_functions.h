#ifndef UTILITY_FUNCTION_H
#define UTILITY_FUNCTION_H

// class containing energy spectrums and other utility functions required by main code

#include "TF1.h"

class utility_functions {

private:
	TF1 *fScintillation_function_muon = nullptr;
	TF1 *fScintillation_function_electron = nullptr;
	TF1 *fScintillation_function_alpha = nullptr;

public:

	// constructor
	utility_functions(){};

	// destructor 
	~utility_functions();

	// poisson distribution function
	int poisson(const double mean, const double draw, const double eng) const;

	// scintillation function
	static double scintillation_function( const double *t, const double *par);
	void initalise_scintillation_functions(const double t_singlet, const double t_triplet, const double scint_time_window);
	double get_scintillation_time(const int &particle_type);

};

#endif