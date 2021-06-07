// main driver for scintillation toy mc simulation code

#include<string>
#include<iostream>
#include<fstream>
#include<chrono>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TRandom.h"
#include "TVector3.h"

#include "data_output.h"
#include "semi_analytic_hits.h"
#include "time_parameterisation.h"
#include "utility_functions.h"
#include "radiological_parameters.h"

// include parameter file
#include "simulation_parameters.h"
using namespace std;

int main() {

	gRandom->SetSeed(0);

	// ------- Initialise output class ----------
	data_output output_file(parameters::output_file_name, parameters::include_timings, parameters::include_reflected);

	// -------- Initialise semi-analytic hits class ---------
	semi_analytic_hits hits_model;

	// -------- Initialise timing parametrisation class ---------
	time_parameterisation times_model(parameters::timing_discretisation_step_size);

	// -------- Initialise utility/energy spectrum class ---------
	utility_functions utility;
	utility.initalise_scintillation_functions_argon(parameters::t_singlet, parameters::t_triplet, parameters::singlet_fraction_electron, parameters::triplet_fraction_electron,
                                                  parameters::singlet_fraction_alpha, parameters::triplet_fraction_alpha, parameters::scint_time_window);
  utility.initalise_scintillation_functions_xenon(parameters::t_singlet_Xe, parameters::t_triplet_Xe, parameters::singlet_fraction_Xe, parameters::triplet_fraction_Xe,
                                                  parameters::scint_time_window);

  // ------- Read photon detector positions and types --------
  std::vector<std::vector<int>> opdet_type;
  std::vector<std::vector<double>> opdet_position;

  std::cout << "Loading Photon Detector positions..." << std::endl;
  std::ifstream detector_positions_file;
  detector_positions_file.open("dunevd_4x7_base_mapping.txt");
  if(detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
  else {std::cout << "File not found." << std::endl; exit(1);}
  while(!detector_positions_file.eof()) {
      int num_opdet, type_opdet, orientation_opdet; double x_opdet, y_opdet, z_opdet;
      if(detector_positions_file >> num_opdet >> x_opdet >> y_opdet >> z_opdet >> type_opdet >> orientation_opdet) {
          std::vector<int> type({num_opdet, type_opdet, orientation_opdet});
          std::vector<double> position({x_opdet, y_opdet, z_opdet});         
          opdet_type.push_back(type);
          opdet_position.push_back(position);
      }
      else{ break; }
  }
  detector_positions_file.close();
  int number_opdets = opdet_type.size();
  std::cout << "Positions Loaded: " << number_opdets << " optical detectors." << std::endl << std::endl;


        //---Alpha-gamma parameters
        double a_gamma_distance;
        double gamma_length = 12.0; //Taken as gamma shower length in LAr, 12 cm
        double added_x;
        double added_y;
        double added_z;

	TFile *f_alpha;
	TH1D *alpha_gamma;
	if(parameters::gen_alpha_gamma == true) {
        //---Reading in positions of alpha decays from .root
        f_alpha = new TFile("../spectra/real_alpha_positions_363.root");
        alpha_gamma = (TH1D*)f_alpha->Get("histo");
	}

        int max_events;
        int scint_yield;
        string particle;


	/////////////-------------Setting beta and gamma spectra---------------////////////////////
        double Q_beta_endpoint = 0;    
        if(parameters::gen_argon == true) {Q_beta_endpoint = radiological::Q_Ar;}
        if(parameters::gen_Co60B == true) {Q_beta_endpoint = radiological::Q_Co60B;}
        if(parameters::gen_Ar42 == true) {Q_beta_endpoint = radiological::Q_Ar42;}
        if(parameters::gen_K42 == true) {Q_beta_endpoint = radiological::Q_K42;}
        if(parameters::gen_40KB == true) {Q_beta_endpoint = radiological::Q_40KB;}
        if(parameters::gen_Kr85B1 == true) {Q_beta_endpoint = radiological::Q_Kr85B1;}
        if(parameters::gen_Kr85B2 == true) {Q_beta_endpoint = radiological::Q_Kr85B2;}
        if(parameters::gen_Pb214 == true) {Q_beta_endpoint = radiological::Q_Pb214;}
        if(parameters::gen_Bi214 == true) {Q_beta_endpoint = radiological::Q_Bi214;}
	//---Gammas---// (Could add some smearing to these spectra, but don't currently)
        if(parameters::gen_Co60G1 == true) {parameters::fixed_energy = true; parameters::fixedE = 1.173;}
        else if(parameters::gen_Co60G2 == true) {parameters::fixed_energy = true; parameters::fixedE = 1.332;}
        else if(parameters::gen_40KG == true) {parameters::fixed_energy = true; parameters::fixedE = 1.46;}
        else if(parameters::gen_Kr85G1 == true) {parameters::fixed_energy = true; parameters::fixedE = 0.151;}
        else if(parameters::gen_Kr85G2 == true) {parameters::fixed_energy = true; parameters::fixedE = 0.305;}


        TF1 *fSpectrum = new TF1("fSpectrum",utility_functions::SpectrumFunction,0,Q_beta_endpoint,1);
        TF1 *flandau_sn = new TF1("flandau_sn",utility_functions::fsn, 0, 50, 1);
        TF1 *flandau_so = new TF1("flandau_so",utility_functions::fso, 0, 16.56, 1); 
        TF1 *flandau_hep = new TF1("flandau_hep",utility_functions::fhep, 0, 18.78, 1); 

        flandau_sn->SetParameter(0, radiological::Eav);
        flandau_so->SetParameter(0, radiological::Eav);
        flandau_hep->SetParameter(0, radiological::Eav);

        TRandom3 *fGauss = new TRandom3();

        //---GAMMA backgrounds
        if(parameters::fixed_energy == true) {
           max_events = parameters::max_events_FE;
        if(parameters::gen_Co60G1 == true) {max_events = parameters::max_events_Co60G1;}
        if(parameters::gen_Co60G2 == true) {max_events = parameters::max_events_Co60G2;}
        if(parameters::gen_40KG == true) {max_events = parameters::max_events_40KG;}
        if(parameters::gen_Kr85G1 == true) {max_events = parameters::max_events_Kr85G1;}
        if(parameters::gen_Kr85G2 == true) {max_events = parameters::max_events_Kr85G2;}

        scint_yield = parameters::scintillation_yield;
        particle = "electron";
	std::cout << std::endl << "Generating " << max_events << " events, of fixed energy: " << parameters::fixedE << " MeV." << std::endl;
        }
	//---Other backgrounds
        if(parameters::gen_argon == true) {
           fSpectrum->SetParameter(0, radiological::Q_Ar);
           max_events = parameters::max_events_Ar;
           scint_yield = parameters::scintillation_yield;
           particle = "electron"; //This isn't used to generate anything, just a label which can be printed
           cout << endl << "Generating " << max_events << ", Ar 39 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Co60B == true) {
           fSpectrum->SetParameter(0, radiological::Q_Co60B);
           max_events = parameters::max_events_Co60B;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Co60 beta decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Ar42 == true) {
           fSpectrum->SetParameter(0, radiological::Q_Ar42);
           max_events = parameters::max_events_Ar42;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Ar42 beta decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_40KB == true) {
           fSpectrum->SetParameter(0, radiological::Q_40KB);
           max_events = parameters::max_events_40KB;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", 40K beta decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Kr85B1 == true) {
           fSpectrum->SetParameter(0, radiological::Q_Kr85B1);
           max_events = parameters::max_events_Kr85B1;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Kr85B1 beta decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Kr85B2 == true) {
           fSpectrum->SetParameter(0, radiological::Q_Kr85B2);
           max_events = parameters::max_events_Kr85B2;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << endl << "Generating " << max_events << ", Kr85B2 beta decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::supernova == true){
           max_events = parameters::max_events_SN;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", supernova events.\n";
        }
        if(parameters::solar == true){
           max_events = parameters::max_events_SO;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", solar events.\n";
        }
        if(parameters::gen_Pb214 == true){
           fSpectrum->SetParameter(0, radiological::Q_Pb214);
           max_events = parameters::max_events_Pb214;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", Pb214 events.\n";
        }
        if(parameters::gen_Bi214 == true){
           fSpectrum->SetParameter(0, radiological::Q_Bi214);
           max_events = parameters::max_events_Bi214;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", Bi214 events.\n";
        }
        if(parameters::gen_K42 == true){
           fSpectrum->SetParameter(0, radiological::Q_K42);
           max_events = parameters::max_events_K42;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", K42 events.\n";
        }
        if(parameters::gen_hep == true){ //hep solar neutrinos
           max_events = parameters::max_events_hep;
           scint_yield = parameters::scintillation_yield;
           particle = "electron";
           cout << "\nGenerating " << max_events << ", hep events.\n";
        }
        if(parameters::gen_Po210 == true){
           max_events = parameters::max_events_Po210;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Po210 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Rn222 == true){
           max_events = parameters::max_events_Rn222;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Rn222 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Po218 == true){
           max_events = parameters::max_events_Po218;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Po218 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_Po214 == true){
           max_events = parameters::max_events_Po214;
           scint_yield = parameters::scint_yield_alpha;
           particle = "alpha";
           cout << "\nGenerating " << max_events << ", Po214 decays in time window: " << parameters::time_window << " seconds." << endl;
        }
        if(parameters::gen_alpha_gamma == true){
           max_events = parameters::max_events_alpha_gamma;
           scint_yield = parameters::scint_yield_alpha;
           cout << "\nGenerating " << max_events << ", Alpha-gamma decays in time window: " << parameters::time_window << " seconds." << endl;
        }


    // ----------- Create Events ------------//
    //generate event positions and energies, storing information in output file	
    std::vector<double> energy_list; 
    energy_list.reserve(max_events);
    std::vector<std::vector<double>> position_list(max_events, std::vector<double>(3,0.0));

    std::cout << "Generating events..." << std::endl;

    for (int event = 0; event < max_events; event++){//Start of event loop

	//--Printing Completion %
        if ( (event != 0) && (max_events >= 10) &&  (event % (max_events/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", event / (max_events/10));
        }

        //--Energy of event
        double energy;
        if(parameters::fixed_energy == true) {energy = parameters::fixedE;}
        if(parameters::gen_argon == true) {energy = fSpectrum->GetRandom();}    
        if(parameters::gen_Ar42 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_K42 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Bi214 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Pb214 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Co60B == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_40KB == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Kr85B1 == true) {energy = fSpectrum->GetRandom();}
        if(parameters::gen_Kr85B2 == true) {energy = fSpectrum->GetRandom();}   
        if(parameters::supernova == true) {energy = flandau_sn->GetRandom();}   
        if(parameters::solar == true) {energy = flandau_so->GetRandom();}   
        if(parameters::gen_hep == true) {energy = flandau_hep->GetRandom();}
        if(parameters::gen_Po210 == true) {energy = fGauss->Gaus(radiological::Q_Po210, 0.05);}
        if(parameters::gen_Rn222 == true) {energy = fGauss->Gaus(radiological::Q_Rn222, 0.05);}
        if(parameters::gen_Po218 == true) {energy = fGauss->Gaus(radiological::Q_Po218, 0.05);}
        if(parameters::gen_Po214 == true) {energy = fGauss->Gaus(radiological::Q_Po214, 0.05);}
	//--Alpha-gamma energy
	if(parameters::gen_alpha_gamma == true) {
	    if(event % 2 ==0) { //Alpha event
                energy = fGauss->Gaus(radiological::Q_Rn222, 0.05);
	    }
            else { //Gamma event
		energy = fGauss->Gaus(15.0, 2.9);
	    }
	}
        energy_list.push_back(energy);

        //--Position of event
        if(parameters::fixed_energy == true || parameters::gen_argon == true || parameters::gen_Ar42 == true || parameters::gen_K42 == true || parameters::gen_Bi214 == true || parameters::gen_Pb214 == true || parameters::gen_Kr85B1 == true || parameters::gen_Kr85B2 == true || parameters::supernova == true || parameters::solar == true || parameters::gen_hep == true || parameters::gen_Kr85G1 == true || parameters::gen_Kr85G2 == true){
           position_list[event][0] = gRandom->Uniform(parameters::x_position_range[0],parameters::x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(parameters::y_position_range[0],parameters::y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(parameters::z_position_range[0],parameters::z_position_range[1]);
        }
        else if(parameters::gen_40KB == true || parameters::gen_40KG == true){
           position_list[event][0] = gRandom->Uniform(radiological::K_x_position_range[0],radiological::K_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(radiological::K_y_position_range[0],radiological::K_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(radiological::K_z_position_range[0],radiological::K_z_position_range[1]);
        }
        else if(parameters::gen_Co60B == true){
           position_list[event][0] = gRandom->Uniform(radiological::Co_x_position_range[0],radiological::Co_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(radiological::Co_y_position_range[0],radiological::Co_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(radiological::Co_z_position_range[0],radiological::Co_z_position_range[1]);
        }
        else if(parameters::gen_Po210 == true){
           position_list[event][0] = gRandom->Uniform(radiological::Po_x_position_range[0],radiological::Po_x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(radiological::Po_y_position_range[0],radiological::Po_y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(radiological::Po_z_position_range[0],radiological::Po_z_position_range[1]);
        }
        else {
           position_list[event][0] = gRandom->Uniform(parameters::x_position_range[0],parameters::x_position_range[1]);
           position_list[event][1] = gRandom->Uniform(parameters::y_position_range[0],parameters::y_position_range[1]);
           position_list[event][2] = gRandom->Uniform(parameters::z_position_range[0],parameters::z_position_range[1]);
        }

        //--Position of alpha-gammas from Rn222 chain
        if(parameters::gen_alpha_gamma == true) {
	    if(event % 2 ==0) { //Alpha event
                position_list[event][0] = alpha_gamma->GetRandom();
                position_list[event][1] = gRandom->Uniform(parameters::y_position_range[0],parameters::y_position_range[1]);
                position_list[event][2] = gRandom->Uniform(parameters::z_position_range[0],parameters::z_position_range[1]);
	    }
	    else { //Gamma event
		a_gamma_distance = gRandom->Exp(gamma_length);
		TRandom obj;
		obj.SetSeed(0);
		obj.Sphere(added_x, added_y, added_z, a_gamma_distance);
		position_list[event][0] = position_list[event-1][0] + added_x;
                position_list[event][1] = position_list[event-1][1] + added_y;
                position_list[event][2] = position_list[event-1][2] + added_z;
            }
	}
                //These statements are to stop gammas exiting through the edges of the TPC
                if (position_list[event][0] > parameters::max_x) {
                    position_list[event][0] = parameters::max_x;
                }
                if (position_list[event][0] < parameters::min_x) {
                    position_list[event][0] = parameters::min_x;
                }
                if (position_list[event][1] > parameters::max_y) {
                    position_list[event][1] = parameters::max_y;
                }
                if (position_list[event][1] < parameters::min_y) {
                    position_list[event][1] = parameters::min_y;
                }
                if (position_list[event][2] > parameters::max_z) {
                    position_list[event][2] = parameters::max_z;
                }
                if (position_list[event][2] < parameters::min_z) {
                    position_list[event][2] = parameters::min_z;
                }


	//---------SET THESE FOR FIXED POSITIONS----------//
        //position_list[event][0] = gRandom->Uniform(310.0,parameters::x_position_range[1]);
        //position_list[event][0] = 5.0; //Fixed x position
        //position_list[event][1] = 0.0; //Fixed y position
        //position_list[event][2] = 200.0; //Fixed z position

        // add event properties to output file
        output_file.add_event(event, energy_list[event], position_list[event]);

    }//End of event loop

	std::cout << "Event generation complete." << std::endl << std::endl;
	
	// --------- Calculate hits and times ----------
	std::cout << "Determining number of photon hits..." << std::endl;

	//--loop each event in events list
	for(int event = 0; event < max_events; event++) {

	  //--output completion %
    if ( (event != 0) && (max_events >= 10) &&  (event % (max_events/10) == 0) ) {
        std::cout << Form("%i0%% Completed...\n", event / (max_events/10));
    }

    // particle type 
    bool isAlpha = false;
    if (parameters::gen_alpha_gamma == true && event % 2 == 0) isAlpha = true;
    else if (parameters::particle_type == 1) isAlpha = true;

    // number of photons produced
	  int number_photons;
	  if(isAlpha) { 
		  number_photons = utility.poisson(static_cast<double>(parameters::scint_yield_alpha) * energy_list.at(event), gRandom->Uniform(1.), energy_list.at(event));
	  }
    else number_photons = utility.poisson(static_cast<double>(parameters::scintillation_yield) * energy_list.at(event), gRandom->Uniform(1.), energy_list.at(event));

    // singlet/triplet fraction
    double singlet_fraction;
    double triplet_fraction;
    if (isAlpha) {
      singlet_fraction = parameters::singlet_fraction_alpha;
      triplet_fraction = parameters::triplet_fraction_alpha;
    }
    else {
      singlet_fraction = parameters::singlet_fraction_electron;
      triplet_fraction = parameters::triplet_fraction_electron;      
    }

    // determine number of photons detected and their timing:
    // get scintillation position
    TVector3 ScintPoint(position_list[event][0],position_list[event][1],position_list[event][2]);
    
    // detector global QE
    double globalQE_VUV = parameters::quantum_efficiency * parameters::wireplane_factor * parameters::vuv_transmission * (parameters::opdet_fraction_both + parameters::opdet_fraction_vuv_only);
    double globalQE_VIS = parameters::quantum_efficiency * parameters::wireplane_factor * parameters::cathode_tpb_frac * parameters::vis_transmission * (parameters::opdet_fraction_both +  parameters::opdet_fraction_visible_only);
    
    // loop over each optical channel
    for(int op_channel = 0; op_channel < number_opdets; op_channel++) {

      // get optical detector type - rectangular or disk aperture
      int op_channel_type = opdet_type[op_channel][1];

      // get optical detector orientation - cathode (x) == 0, lateral (y) == 1
      int op_channel_orientation = opdet_type[op_channel][2];

      // get detection channel coordinates (in cm)
      TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);

      // determine number of hits on optical channel via semi-analytic model:  

      // VUV
      // calculate detected photons
      int num_VUV_Ar = 0;
      int num_VUV_Xe = 0;
      if (parameters::simulate_xenon == false) { // argon only case            
        // incident photons
        int num_VUV_geo = hits_model.VUVHits(number_photons, ScintPoint, OpDetPoint, op_channel_type, op_channel_orientation, 0);       // calculate hits       
        // apply additional factors QE etc.
        double detectorQE_VUV = globalQE_VUV;
        if (op_channel_orientation == 1) detectorQE_VUV = globalQE_VUV * parameters::field_cage_transparency; //laterals                
        for(int i = 0; i < num_VUV_geo; i++) if (gRandom->Uniform(1.) <= detectorQE_VUV) num_VUV_Ar++;
      }
      if (parameters::simulate_xenon == true) { // xenon doped case         
        // split into prompt and late light
        int number_photons_Ar = std::round(number_photons*singlet_fraction);
        int number_photons_Xe = std::round(number_photons*triplet_fraction);

        // incident photons
        int num_VUV_geo_Ar = hits_model.VUVHits(number_photons_Ar, ScintPoint, OpDetPoint, op_channel_type, op_channel_orientation, 0);       // prompt light as argon
        int num_VUV_geo_Xe = hits_model.VUVHits(number_photons_Xe, ScintPoint, OpDetPoint, op_channel_type, op_channel_orientation, 1);       // late light as xenon       
        // apply additional factors QE etc.
        double detectorQE_VUV = globalQE_VUV;
        if (op_channel_orientation == 1) detectorQE_VUV = globalQE_VUV * parameters::field_cage_transparency; //laterals            
        for(int i = 0; i < num_VUV_geo_Ar; i++) if (gRandom->Uniform(1.) <= detectorQE_VUV) num_VUV_Ar++;
        for(int i = 0; i < num_VUV_geo_Xe; i++) if (gRandom->Uniform(1.) <= detectorQE_VUV) num_VUV_Xe++;
      }  

      // Visible (foils)
      // calculate detected photons
      int num_VIS_Ar = 0;
      int num_VIS_Xe = 0;
      if (parameters::include_reflected) {
        if (parameters::simulate_xenon == false) { // argon only case  
          // incident photons
          int num_VIS_geo = hits_model.VisHits(number_photons, ScintPoint, OpDetPoint, op_channel_type, 0);     // calculate hits       
          // apply additional factors QE etc.
          for(int j = 0; j < num_VIS_geo; j++) if (gRandom->Uniform(1.) <= globalQE_VIS) num_VIS_Ar++;
        }
        if (parameters::simulate_xenon == true) { // xenon doped case         
          // split into prompt and late light
          int number_photons_Ar = std::round(number_photons*singlet_fraction);
          int number_photons_Xe = std::round(number_photons*triplet_fraction);

          // incident photons
          int num_VIS_geo_Ar = hits_model.VisHits(number_photons_Ar, ScintPoint, OpDetPoint, op_channel_type, 0);     // prompt light as argon
          int num_VIS_geo_Xe = hits_model.VisHits(number_photons_Xe, ScintPoint, OpDetPoint, op_channel_type, 1);     // late light as xenon      
          // apply additional factors QE etc.
          for(int j = 0; j < num_VIS_geo_Ar; j++) if (gRandom->Uniform(1.) <= globalQE_VIS) num_VIS_Ar++;
          for(int j = 0; j < num_VIS_geo_Xe; j++) if (gRandom->Uniform(1.) <= globalQE_VIS) num_VIS_Xe++;
        }
      }
     
      // if no photons from this event for this optical channel, go to the next channel.
      int num_VUV = num_VUV_Ar + num_VUV_Xe;
      int num_VIS = num_VIS_Ar + num_VIS_Xe;
      if(num_VUV+num_VIS == 0) { continue; } // forces the next iteration

      // calculate timings
      std::vector<double> total_time_vuv; total_time_vuv.reserve(num_VUV);
      std::vector<double> total_time_vis; total_time_vis.reserve(num_VIS);
      if (parameters::include_timings){
        // split into Ar/Xe
        std::vector<double> total_time_vuv_Ar; total_time_vuv_Ar.reserve(num_VUV_Ar);
        std::vector<double> total_time_vuv_Xe; total_time_vuv_Xe.reserve(num_VUV_Xe);
        std::vector<double> total_time_vis_Ar; total_time_vis_Ar.reserve(num_VIS_Ar);
        std::vector<double> total_time_vis_Xe; total_time_vis_Xe.reserve(num_VIS_Xe);

        // VUV, Ar                  
        if(num_VUV_Ar > 0) {
          
          // transport times
          double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();
          double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance_to_pmt;
          double theta = acos(cosine)*180./3.14159;
          int angle_bin = theta/45;       // 45 deg bins    

          std::vector<double> transport_time_vuv_Ar = times_model.getVUVTime(distance_to_pmt, angle_bin, num_VUV);

          for(auto& x: transport_time_vuv_Ar) {
      			// emission time
            double emission_time;
            if (parameters::simulate_xenon == false) {
              if (isAlpha) emission_time = utility.get_scintillation_time_alpha()*1000000.0; // in us
              else emission_time = utility.get_scintillation_time_electron()*1000000.0; // in us
            }
            if (parameters::simulate_xenon == true) { // in this case all remaining argon light is prompt
              emission_time = utility.get_scintillation_time_prompt()*1000000.0; // in us
            }
           
            // total time
            double total_time = (x*0.001 + emission_time + 2.5*0.001); // in microseconds	   
            total_time_vuv_Ar.push_back(total_time);
          }
        }
        // VUV, Xe
        if(num_VUV_Xe > 0) {

          // transport times
          double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();
          std::vector<double> transport_time_vuv_Xe = times_model.getVUVTimeXe(distance_to_pmt, num_VUV);

          for(auto& x: transport_time_vuv_Xe) {
            // emission time
            double emission_time = utility.get_scintillation_time_xenon()*1000000.0; // in us       
           
            // total time
            double total_time = (x*0.001 + emission_time + 2.5*0.001); // in microseconds    
            total_time_vuv_Xe.push_back(total_time);
          }
        } 
        
        // VIS,
        if (num_VIS > 0 && parameters::include_reflected) {
          // Ar
          if (num_VIS_Ar > 0) {
            
            // transport times
            std::vector<double> transport_time_vis_Ar = times_model.getVisTime(ScintPoint, OpDetPoint, num_VIS_Ar);

            for(auto& y: transport_time_vis_Ar) {
              // emission time
              double emission_time;
              if (parameters::simulate_xenon == false) {
                if (isAlpha) emission_time = utility.get_scintillation_time_alpha()*1000000.0; // in us
                else emission_time = utility.get_scintillation_time_electron()*1000000.0; // in us
              }
              if (parameters::simulate_xenon == true) { // in this case all remaining argon light is prompt
                emission_time = utility.get_scintillation_time_prompt()*1000000.0; // in us
              }

              // total time
              double total_time = (y*0.001 + emission_time + 2.5*0.001); // in microseconds
              total_time_vis_Ar.push_back(total_time);
            }
          }
          // Xe
          if(num_VIS_Xe > 0) {

            // transport times
            std::vector<double> transport_time_vis_Xe = times_model.getVisTimeXe(ScintPoint, OpDetPoint, num_VIS_Xe);

            for(auto& x: transport_time_vis_Xe) {
              // emission time
              double emission_time = utility.get_scintillation_time_xenon()*1000000.0; // in us       
             
              // total time
              double total_time = (x*0.001 + emission_time + 2.5*0.001); // in microseconds    
              total_time_vis_Xe.push_back(total_time);
            }
          }          
        }

        // combine timings into single vectors for Direct and Reflected light
        total_time_vuv = total_time_vuv_Ar; total_time_vuv.insert( total_time_vuv.end(), total_time_vuv_Xe.begin(), total_time_vuv_Xe.end() );
        total_time_vis = total_time_vis_Ar; total_time_vis.insert( total_time_vis.end(), total_time_vis_Xe.begin(), total_time_vis_Xe.end() );    

      } // end timings block

      // fill data trees for each photon detected
      if (parameters::include_timings) output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint, total_time_vuv, total_time_vis);
      else output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint);

    } // end of optical channel loop
	
	} // end of event loop

  // close alpha file if opened
	if(parameters::gen_alpha_gamma == true) {
	  f_alpha->Close();
	}

	// write output root file
	output_file.write_output_file();

	std::cout << "Program finished." << std::endl;
}
