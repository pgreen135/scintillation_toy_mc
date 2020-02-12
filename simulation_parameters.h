// this file contains parameters used in simulation

namespace parameters {
	
	// output file
	const char *output_file_name = "threading_testing.root";

	// events to generate
	// number
	const int number_events =  10000;
	
	// type

	// energy
	const double energy = 100;	// MeV
	
	// position
	// random position ranges
	const double x_position_range[2] {10,350};	// cm
	const double y_position_range[2] {-600,600};	// cm
	const double z_position_range[2] {400,1000};	// cm


	// semi-analytic simulation options
	// visible light
	const bool include_reflected = true;

	// timings
	const bool include_timings = true;
	const double timing_discretisation_step_size = 1.0;	// cm


	// photon detection system properties
	const double quantum_efficiency = 0.035;				// arapuca QE
	const double mesh_factor = 0.7;	
	const double vuv_tpb_transmission = 0.5;
	const double vis_tpb_transmission = 0.6; 
	const double opdet_tpb_frac = 1.0; 
	const double cathode_tpb_frac = 0.8;		


	// scintillation properties
	const int scintillation_yield = 24000; 		// 24000 photons/MeV at 500 V/m
	const double t_singlet = 0.000000006; 		// 6ns 
	const double t_triplet = 0.0000015; 		// 1.5 us
	const double scint_time_window = 0.00001; 	// 10 us
	const double particle_type = 0;			// ionising particle: 0 - electron, 1 - alpha
		
}
