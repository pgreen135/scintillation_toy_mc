// this file contains parameters used in simulation

#include <string>

namespace parameters {
	
	// output file
	const char *output_file_name = "testing.root";

	// events to generate
	const bool event_from_file = true;
	const std::string event_file_name = "sbnd_michel_electron.txt"; 
	
	// define events to generate if not from file:
	// number
	const int number_events =  500;	
	// energy
	const double energy = 25;	// MeV	
	// position
	// random position ranges
	const double x_position_range[2] {10,190};	// cm
	const double y_position_range[2] {-190,190};	// cm
	const double z_position_range[2] {10,490};	// cm
	// particle type
	const double particle_type = 1;			// ionising particle: 0 = muon, 1 = electron, 2 = alpha


	// semi-analytic simulation options
	// visible light
	const bool include_reflected = true;

	// timings
	const bool include_timings = true;
	const double timing_discretisation_step_size = 1.0;	// cm


	// photon detection system properties
	const double quantum_efficiency = 0.25;	// PMT QE
	const double wireplane_factor = 1.0;	
	const double vuv_transmission = 1.0;	// efficiency of any detector to VUV photons (128 or 174 nm)
	const double vis_transmission = 1.0; 	// efficiency of any detector to Visible photons
	
	const double opdet_fraction_vuv_only = 0.0;
	const double opdet_fraction_visible_only = 0.0; 
	const double opdet_fraction_both = 1.0 - opdet_fraction_vuv_only - opdet_fraction_visible_only;
	
	const double cathode_tpb_frac = 1.0;		


	// scintillation properties
	const int scintillation_yield_muon = 24000; 		// Muon: 24000 photons/MeV at 500 V/m
	const int scintillation_yield_electron = 20000; 	// Electron: 20000 photons/MeV at 500 V/m
	const int scintillation_yield_alpha = 20000; 		// Alpha: 20000 photons/MeV at 500 V/m
	const double t_singlet = 0.000000006; 		// 6ns 
	const double t_triplet = 0.0000015; 		// 1.5 us
	const double scint_time_window = 0.00001; 	// 10 us
	
		
}
