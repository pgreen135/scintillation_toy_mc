// this file contains parameters used in simulation

#include "TRandom3.h"
namespace parameters {
	
	// output file
	const char *output_file_name = "testing.root";


    ///-------------------------------------
    //--------Simulation Settings-----------
    ///-------------------------------------
    // Xenon doping
    const bool simulate_xenon = true;   // enables xenon doping, false = pure argon

    // WLS Reflective foils
    const bool include_reflected = false; // enables WLS relfective foils on the cathode (visible light) -- not applicable for VD

    // Timings
    const bool include_timings = true;   // enables timings (emission and transport)  


    ///-------------------------------------
    //--------WHAT events to generate?------
    ///-------------------------------------
    bool fixed_energy = true; //If gamma activity this will be set to true in the analyze_light.cpp
    double fixedE = 10.0;
    
    const double particle_type = 0;   // ionising particle type: 0 - electron, 1 - alpha

    bool supernova = false;
    bool solar = false;
    bool gen_hep = false;
    bool gen_argon = false; //This is Ar39
    bool gen_Co60B = false;
    bool gen_Co60G1 = false; //fixed energy = 1.173 Mev
    bool gen_Co60G2 = false; //fixed energy = 1.332 Mev
    bool gen_Ar42 = false;
    bool gen_K42 = false;
    bool gen_40KB = false;
    bool gen_40KG = false; //fixed energy = 1.46 MeV
    bool gen_Kr85B1 = false;
    bool gen_Kr85B2 = false;
    bool gen_Kr85G1 = false; //fixed energy = 0.151 Mev
    bool gen_Kr85G2 = false; //fixed energy = 0.305 Mev
    bool gen_Pb214 = false;
	//////-------------------------------///////
	//////--Radon alpha chain daughters--///////
	//////-------------------------------///////
	bool gen_Rn222 = false;
	bool gen_Po218 = false;
	bool gen_Po214 = false;
	bool gen_Po210 = false;
    bool gen_Bi214 = false;
	bool gen_alpha_gamma = false; //specifically for generating an alpha event followed by a rare gamma emission

	//---alpha scintillation properties---//
    const double massAlpha = 3727.3794;     // alpha particle mass - MeV/c^2
	const double activity_Rn = 0.000021;    // Bq/kg
	const double activity_Ar = 1.01;

	//----------------------------------------------
	//---------Dimensions of event region-----------
	//----------------------------------------------
	// region to generate events in
    const double x_position_range[2] {-320,320};	// cm 
	const double y_position_range[2] {-600,600};	// cm
	const double z_position_range[2] {10,2000};   // cm    
    
    //-------------------------------------
    //----Number of Events-----------------
    //-------------------------------------
    // Fixed energy (electron like) events:
    const int max_events_FE = 10000;

    //chrisflynn mass
    const double entire_mass = 0.001396*(x_position_range[1]-x_position_range[0])*(y_position_range[1]-y_position_range[0])*(z_position_range[1]-z_position_range[0]); //density = 1.396 g/cm^3

    //chrisflynn random variables
    const double frame_time = 0.0025; // 2.5 miliseconds; equivalent to the drift time of electrons in DUNE
    const double time_window = 1; //readout window (running time of the simulation in seconds)
    const double time_frames = int(time_window/frame_time); // number of frames registered
    
	// Ar-39 events:
    //const int max_events_Ar = activity_Ar * entire_mass * time_window;//FULL volume for 1 TPC
	const int max_events_Ar = 1000;
    const int Ar_decays_per_sec = activity_Ar* entire_mass; // decay rate in one TPC
        
	// Radon events:
    const int max_events_Rn = 5.584e-5 * (x_position_range[1]-x_position_range[0])*(y_position_range[1]-y_position_range[0])*(z_position_range[1]-z_position_range[0]) * time_window * 0.01; //the 0.01 is included because of design changes to reduce alpha see https://indico.fnal.gov/event/20794/contribution/0/material/slides/0.pdf
    const int max_events_Rn222 = max_events_Rn * 0.25;
	const int max_events_Po218 = max_events_Rn * 0.25; //0.25 comes from 4 alpha decays in chain from Rn222, each contributes 10Bq/kg * 0.01 with new Rn222 condition
    const int max_events_Po214 = max_events_Rn * 0.25;
	//Will want a new variable for Po210 here if including events from Rn222 decay chain 
    const int max_events_Po210 = 5e-6 * (radiological::Po_x_position_range[1]-radiological::Po_x_position_range[0])*(radiological::Po_y_position_range[1]-radiological::Po_y_position_range[0])*(radiological::Po_z_position_range[1]-radiological::Po_z_position_range[0]) * time_window;//Activity from Po210 in the APAs
    const int max_events_Pb214 = max_events_Rn * 0.25;  
    const int max_events_Bi214 = max_events_Rn * 0.25;  
    const double Rn_decays_per_sec = activity_Rn* entire_mass; // decay rate in one TPC
        
	//Supernova events:
    const int max_events_SN = time_frames;
        
	//8B Solar neutrino events:
    const int max_events_SO = 1000;

	// Other radiactive bkgs
    const int max_events_Co60B = 8.2e-5 * (radiological::Co_x_position_range[1]-radiological::Co_x_position_range[0])*(radiological::Co_y_position_range[1]-radiological::Co_y_position_range[0])*(radiological::Co_z_position_range[1]-radiological::Co_z_position_range[0]) * time_window;//For beta decay, near 100% would have beta decay and 200% gamma decay
    const int max_events_Co60G1 = max_events_Co60B;
    const int max_events_Co60G2 = max_events_Co60B;
    const int max_events_Ar42 = 1.283768e-7 * (x_position_range[1]-x_position_range[0])*(y_position_range[1]-y_position_range[0])*(z_position_range[1]-z_position_range[0]) * time_window;//1.41e-3 Bq*cm^-3 from mcc11 simulation
    const int max_events_K42 = max_events_Ar42 * 0.819;
	const int max_events_40KB = 10000000;
    const int max_events_40KG = 2.7195e-3 * (radiological::K_x_position_range[1]-radiological::K_x_position_range[0])*(radiological::K_y_position_range[1]-radiological::K_y_position_range[0])*(radiological::K_z_position_range[1]-radiological::K_z_position_range[0]) * time_window * 0.1072;//10.72% ratio for gamma
    const int max_events_Kr85B1 = 1.6e-4 * (x_position_range[1]-x_position_range[0])*(y_position_range[1]-y_position_range[0])*(z_position_range[1]-z_position_range[0]) * time_window * 0.785;//near 100% of Kr82 would have beta decay
    const int max_events_Kr85B2 = 1.6e-4 * (x_position_range[1]-x_position_range[0])*(y_position_range[1]-y_position_range[0])*(z_position_range[1]-z_position_range[0]) * time_window * 0.14;
    const int max_events_Kr85G1 = 1.6e-4 * (x_position_range[1]-x_position_range[0])*(y_position_range[1]-y_position_range[0])*(z_position_range[1]-z_position_range[0]) * time_window * 0.785 * 0.752;
    const int max_events_Kr85G2 = max_events_Kr85B2;
    const int max_events_hep = 100;
	
	// Alpha-gamma events
    const int max_events_alpha_gamma = 2000; //This is the number of alpha and gammas added together, e.g for 2 events get 1 alpha and 1 gamma

    // TPC boundaries (for alpha-gammas)
    const double max_x = 363.0;
    const double min_x = -363.0;
    const double max_y = 600.0;
    const double min_y = -600.0;
    const double max_z = 2000.0;
    const double min_z = 0.0;

    // timing parametersiation properties
    const double timing_discretisation_step_size = 1.0; // cm
    
	// photon detection system properties
	const double quantum_efficiency = 0.035;	   // arapuca QE
	const double wireplane_factor = 1.0;           // included in parameterisations	
	const double vuv_transmission = 0.5;           // PD window transmission: VUV Ar/Xe scintillation light
	const double vis_transmission = 0.7;           // PD window transmission: Visible light (foils)
    
    // fraction PDs sensitive to each component
    const double opdet_fraction_vuv_only = 0.0;
    const double opdet_fraction_visible_only = 0.0;
    const double opdet_fraction_both = 1.0 - opdet_fraction_vuv_only - opdet_fraction_visible_only;
 
    // cathode foils/TPB coverage
	const double cathode_tpb_frac = 0.8;

    const double field_cage_transparency = 0.7; // transparency of field cage for lateral PDs       

	// scintillation properties
	const int scintillation_yield = 24000; 		// 24000 photons/MeV at 500 V/m
    const int scint_yield_alpha = 16800;        // SY of alpha particles at 500 V/cm - from larsoft  

    // scintillation timing properties
    const double scint_time_window = 0.00001;   // 10 us
    // argon
    // timing
    const double t_singlet = 0.000000006;       // 6ns 
    const double t_triplet = 0.0000015;         // 1.5 us
    // prompt/late ratio
    // electron-like
    const double singlet_fraction_electron = 0.30;
    const double triplet_fraction_electron = 0.70;
    // alpha
    const double singlet_fraction_alpha = 0.75;
    const double triplet_fraction_alpha = 0.25;   
	
    // xenon -- Aprile, Doke https://arxiv.org/abs/0910.4956
    // timing
    const double t_singlet_Xe = 0.0000000042;      // 4.2 ns 
    const double t_triplet_Xe = 0.000000022;       // 22 ns
    // prompt/late ratio
    // minimal dependence on particle type, approximate as single value
    const double singlet_fraction_Xe = 0.30;
    const double triplet_fraction_Xe = 0.70;    
		
}
