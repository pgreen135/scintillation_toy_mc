// main driver for scintillation toy mc simulation code

#include<string>
#include<iostream>
#include<fstream>
#include<chrono>

#include "TRandom.h"
#include "TVector3.h"

#include "data_output.h"
#include "semi_analytic_hits.h"
#include "time_parameterisation.h"
#include "utility_functions.h"


// include parameter file
#include "simulation_parameters.h"

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
	utility.initalise_scintillation_functions(parameters::t_singlet, parameters::t_triplet, parameters::scint_time_window);
    
	// ------- Read photon detector positions and types --------
	std::vector<std::vector<int>> opdet_type;
	std::vector<std::vector<double>> opdet_position;

	std::cout << "Loading Photon Detector positions..." << std::endl;
    std::ifstream detector_positions_file;
    detector_positions_file.open("pmts_sbnd_v01_04_modified.txt");
    if(detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
    else {std::cout << "File not found." << std::endl; exit(1);}
    while(!detector_positions_file.eof()) {
        int num_opdet, type_opdet; double x_opdet, y_opdet, z_opdet;
        if(detector_positions_file >> num_opdet >> x_opdet >> y_opdet >> z_opdet) {
            if (std::abs(x_opdet) == 206.465) type_opdet = 1; // disk PMT
            else type_opdet = 0; // rectangular Arapuca
            std::vector<int> type({num_opdet, type_opdet});
            std::vector<double> position({x_opdet, y_opdet, z_opdet});         
            opdet_type.push_back(type);
            opdet_position.push_back(position);
        }
        else{ break; }
    }
    detector_positions_file.close();
    int number_opdets = opdet_type.size();
    std::cout << "Positions Loaded: " << number_opdets << " optical detectors." << std::endl << std::endl;

	// ----------- Create Events ------------
	// generate event positions and energies, storing information in output file	

	std::vector<double> energy_list; energy_list.reserve(parameters::number_events); // event energy in MeV
    std:std::vector<double> time_list; time_list.reserve(parameters::number_events); // event time in ns
    std::vector<double> particle_type_list; particle_type_list.reserve(parameters::number_events); // event time in ns	
    std::vector<std::vector<double>> position_list(parameters::number_events, std::vector<double>(3,0.0)); 

	std::cout << "Generating events..." << std::endl;
	if (parameters::event_from_file) {
        int event_counter = 0;
        position_list.clear();
        
        // ------- Read events from text file --------
        
        std::cout << "Reading events from file: " << std::endl;
        std::ifstream event_file;
        event_file.open(parameters::event_file_name);
        if(event_file.is_open()) std::cout << "File opened successfully" << std::endl;
        else {std::cout << "File not found." << std::endl; exit(1);}
        while(!event_file.eof()) {
            int number, pdg; double time, x, y, z; double energy; 
            if(event_file >> number >> pdg >> time >> x >> y >> z >> energy) {
                std::cout << number << " " << pdg << " " << time << " " << x << " " << y << " " << z << " " << energy << std::endl;
                // populate lists with event
                // type
                int particle_type;
                if (pdg == 13 || pdg == -13) particle_type = 0; // muon
                else if (pdg == 11 || pdg == -11) particle_type = 1; // electron
                else particle_type = -1; // unknown
                particle_type_list.push_back(particle_type);

                // initial time
                time_list.push_back(time);
                
                // position
                std::vector<double> position = {x, y, z};
                position_list.push_back(position);
                
                // energy
                // HACKED FOR PARTICULAR EVENT, NEED TO GENERALISE
                // energy per step of track, constant dEdx
                double energy_step;
                if (pdg == 13 || pdg == -13) energy_step = energy/47; // muon
                else if (pdg == 11 || pdg == -11) energy_step = energy/17; // muon
                else energy_step = 0;
                energy_list.push_back(energy_step);
                             
                // add event properties to output file
                output_file.add_event(event_counter, energy_step, time, particle_type, position);
                event_counter++;
               
            }
            else{ break; }
        }
        event_file.close();
        int number_events = energy_list.size();
        std::cout << "Loaded: " << number_events << " events from file." << std::endl << std::endl;
    }
    else {
        for (int event = 0; event < parameters::number_events; event++){

    		// output completion %
            if ( (event != 0) && (parameters::number_events >= 10) &&  (event % (parameters::number_events/10) == 0) ) {
                std::cout << Form("%i0%% Completed...\n", event / (parameters::number_events/10));
            }

            // determine energy of event
            energy_list.push_back(parameters::energy);

            // initial time
            time_list.push_back(0);     // to be implemented

            // particle type;
            particle_type_list.push_back(parameters::particle_type);    // 0 = muon, 1 = electron, 2 = alpha (not fully implemented)

            // determine position of event
            // choose random position in 3D range
            position_list[event][0] = gRandom->Uniform(parameters::x_position_range[0],parameters::x_position_range[1]);
            position_list[event][1] = gRandom->Uniform(parameters::y_position_range[0],parameters::y_position_range[1]);
            position_list[event][2] = gRandom->Uniform(parameters::z_position_range[0],parameters::z_position_range[1]);

            // add event properties to output file
            output_file.add_event(event, energy_list[event], time_list[event], particle_type_list[event], position_list[event]);
    	}
    }
	std::cout << "Event generation complete." << std::endl << std::endl;
	
	// --------- Calculate hits and times ----------
	std::cout << "Determining number of photon hits..." << std::endl;

	// loop each event in events list
    int number_events = 0;
    if (parameters::event_from_file) number_events = energy_list.size();
    else number_events = parameters::number_events;
	for(int event = 0; event < number_events; event++) {

		// output completion %
        if ( (event != 0) && (number_events >= 10) &&  (event % (number_events/10) == 0) ) {
            std::cout << Form("%i0%% Completed...\n", event / (number_events/10));
        }

        // calculate total scintillation yield from the event
        int number_photons = 0;
        // muon
        if (particle_type_list[event] == 0) number_photons = utility.poisson(static_cast<double>(parameters::scintillation_yield_muon) * energy_list.at(event), gRandom->Uniform(1.), energy_list.at(event));
        // electron
        if (particle_type_list[event] == 1) number_photons = utility.poisson(static_cast<double>(parameters::scintillation_yield_electron) * energy_list.at(event), gRandom->Uniform(1.), energy_list.at(event));
        // alpha
        if (particle_type_list[event] == 2) number_photons = utility.poisson(static_cast<double>(parameters::scintillation_yield_alpha) * energy_list.at(event), gRandom->Uniform(1.), energy_list.at(event));
        
        // loop over each optical channel
        for(int op_channel = 0; op_channel < number_opdets; op_channel++) { 

        	// get optical detector type - rectangular or disk aperture
        	int op_channel_type = opdet_type[op_channel][1];

        	// get scintillation point and detection channel coordinates (in cm)
            TVector3 ScintPoint(position_list[event][0],position_list[event][1],position_list[event][2]);
            TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);

            if (OpDetPoint[0] < -50) continue; // 0 photons for opposite TPC

            // determine number of hits on optical channel via semi-analytic model:
            // VUV
            int num_VUV = 0;
            // incident photons
            int num_VUV_geo = hits_model.VUVHits(number_photons, ScintPoint, OpDetPoint, op_channel_type);       // calculate hits       
            // apply additional factors QE etc.            
            for(int i = 0; i < num_VUV_geo; i++) {
                if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::wireplane_factor * parameters::vuv_transmission * (parameters::opdet_fraction_both + parameters::opdet_fraction_vuv_only)) num_VUV++;   
            }
            // Visible
            int num_VIS = 0;
            if (parameters::include_reflected) {
            	// incident photons
            	int num_VIS_geo = 0;                
                num_VIS_geo = hits_model.VisHits(number_photons, ScintPoint, OpDetPoint, op_channel_type);  // calculate hits with hotspot model        
                // apply additional factors QE etc.
            	for(int j = 0; j < num_VIS_geo; j++) {
                	if (gRandom->Uniform(1.) <= parameters::quantum_efficiency * parameters::wireplane_factor * parameters::cathode_tpb_frac * parameters::vis_transmission * (parameters::opdet_fraction_both +  parameters::opdet_fraction_visible_only)) num_VIS++;
                }
            }

            // if no photons from this event for this optical channel, go to the next channel.
            if(num_VUV+num_VIS == 0) { continue; } // forces the next iteration

            // calculate timings
            std::vector<double> total_time_vuv; total_time_vuv.reserve(num_VUV);
            std::vector<double> total_time_vis; total_time_vis.reserve(num_VIS);
            if (parameters::include_timings){
            	// VUV            	
            	if(num_VUV > 0) {
            		// transport times
            		double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();
            		double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance_to_pmt;
  					double theta = acos(cosine)*180./3.14159;
  					int angle_bin = theta/45;	// 45 deg bins	  
            		
            		std::vector<double> transport_time_vuv = times_model.getVUVTime(distance_to_pmt, angle_bin, num_VUV);
            		
            		// total times
            		for(auto& x: transport_time_vuv) {
            			double total_time = (time_list[event] + x + utility.get_scintillation_time(particle_type_list[event])*1e9); // in nanoseconds
            			total_time_vuv.push_back(total_time);
            		}
            	}
            	// VIS
            	if (num_VIS > 0 && parameters::include_reflected) {
            		// transport times
            		std::vector<double> transport_time_vis = times_model.getVisTime(ScintPoint, OpDetPoint, num_VIS);
            		// total times
            		for(auto& y: transport_time_vis) {
            			double total_time = (time_list[event] + y + utility.get_scintillation_time(particle_type_list[event])*1e9 ); // in nanoseconds
            			total_time_vis.push_back(total_time);
            		}
            	}
            }

            // fill data trees for each photon detected
            if (parameters::include_timings) output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint, total_time_vuv, total_time_vis);
            else output_file.add_data(event, op_channel, num_VUV, num_VIS, ScintPoint);

        } // end of optical channel loop

	} // end of event loop

	// write output root file
	output_file.write_output_file();

	std::cout << "Program finished." << std::endl;
}