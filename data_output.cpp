#include "data_output.h"

#include <iostream>

#include "TString.h"

// constructor
data_output::data_output(const char* output_file_name, const bool include_timings, const bool include_reflected) : include_reflected{include_reflected} {

	// create file
	output_file = new TFile(output_file_name, "RECREATE", "Output File");
	if (output_file->IsOpen()) {
		std::cout << "Output file created successfully." << std::endl << std::endl;
	}
	else {
		std::cout << "Output file could not be opened, check file name is valid." << std::endl << std::endl; exit(1);
	} 

	// create trees
	event_tree = new TTree("event_tree", "event tree");
	data_tree = new TTree("data_tree", "data tree");
	data_tree_vuv = new TTree("data_tree_vuv", "data tree_vuv");
	if (include_reflected) data_tree_vis = new TTree("data_tree_vis", "data tree_vis");	

	// create branches
	// event tree
    event_tree->Branch("event_no", &event_no, "event_no/I");
    event_tree->Branch("event_x_pos", &event_x_pos, "event_x_pos/D");
    event_tree->Branch("event_y_pos", &event_y_pos, "event_y_pos/D");
    event_tree->Branch("event_z_pos", &event_z_pos, "event_z_pos/D");
    event_tree->Branch("event_E", &event_E, "event_E/D");
    event_tree->Branch("event_time", &event_time, "event_time/D");
    event_tree->Branch("event_particle_type", &event_particle_type, "event_particle_type/I");
	// data tree containing all light
	if (include_timings) { data_tree->Branch("data_time", &data_time, "data_time/D"); }
    data_tree->Branch("data_pmt", &data_pmt, "data_pmt/I");
    data_tree->Branch("data_x_pos", &data_x_pos, "data_x_pos/D");
    data_tree->Branch("data_y_pos", &data_y_pos, "data_y_pos/D");
    data_tree->Branch("data_z_pos", &data_z_pos, "data_z_pos/D"); 
    data_tree->Branch("data_event", &data_event, "data_event/I");
    // data tree containing vuv (direct) light only
    if (include_timings) { data_tree_vuv->Branch("data_time_vuv", &data_time_vuv, "data_time_vuv/D"); }
    data_tree_vuv->Branch("data_pmt_vuv", &data_pmt_vuv, "data_pmt_vuv/I");
    data_tree_vuv->Branch("data_x_pos_vuv", &data_x_pos_vuv, "data_x_pos_vuv/D"); 
    data_tree_vuv->Branch("data_y_pos_vuv", &data_y_pos_vuv, "data_y_pos_vuv/D");
    data_tree_vuv->Branch("data_z_pos_vuv", &data_z_pos_vuv, "data_z_pos_vuv/D");
    data_tree_vuv->Branch("data_event_vuv", &data_event_vuv, "data_event_vuv/I");
    // data tree containing visible (reflected) light only
    if (include_reflected){
	    if (include_timings) { data_tree_vis->Branch("data_time_vis", &data_time_vis, "data_time_vis/D"); }
	    data_tree_vis->Branch("data_pmt_vis", &data_pmt_vis, "data_pmt_vis/I");
	    data_tree_vis->Branch("data_x_pos_vis", &data_x_pos_vis, "data_x_pos_vis/D");
	    data_tree_vis->Branch("data_y_pos_vis", &data_y_pos_vis, "data_y_pos_vis/D");
	    data_tree_vis->Branch("data_z_pos_vis", &data_z_pos_vis, "data_z_pos_vis/D");
	    data_tree_vis->Branch("data_event_vis", &data_event_vis, "data_event_vis/I");
	}
}

// destructor
data_output::~data_output(){
	// deleting output_file also deletes all trees properly 
	delete output_file;
}

// function that add event information to output file branch
void data_output::add_event(const int &event_number, const double &event_energy, const double &time, const double &particle_type, const std::vector<double> &event_position) {
	event_no = event_number;
	event_E = event_energy; 
	event_time = time;
	event_particle_type = particle_type;
	event_x_pos = event_position[0];
	event_y_pos = event_position[1];
	event_z_pos = event_position[2];
	event_tree->Fill();
}

void data_output::add_data(const int &event_number, const int &optical_channel, const int &num_VUV, const int &num_VIS, const TVector3 &ScintPoint) {
	data_x_pos = ScintPoint[0];
    data_x_pos_vuv = ScintPoint[0];
    data_x_pos_vis = ScintPoint[0];

    data_y_pos = ScintPoint[1];
    data_y_pos_vuv = ScintPoint[1];
    data_y_pos_vis = ScintPoint[1];

    data_z_pos = ScintPoint[2];
    data_z_pos_vuv = ScintPoint[2];
    data_z_pos_vis = ScintPoint[2];

    data_event = event_number;
    data_event_vuv = event_number;
    data_event_vis = event_number;

	data_pmt = optical_channel;
	data_pmt_vuv = optical_channel;
	data_pmt_vis = optical_channel;

	// add entry for each photon
	for (int i = 0; i < num_VUV; i++) {
		data_tree_vuv->Fill();
        data_tree->Fill(); 
	}
	if (include_reflected) {
		for (int i = 0; i < num_VIS; i++) {
			data_tree_vis->Fill();
	        data_tree->Fill(); 
		}
	}
}

void data_output::add_data(const int &event_number, const int &optical_channel, const int &num_VUV, const int &num_VIS, const TVector3 &ScintPoint, const std::vector<double> &times_vuv, const std::vector<double> &times_vis) {
	data_x_pos = ScintPoint[0];
    data_x_pos_vuv = ScintPoint[0];
    data_x_pos_vis = ScintPoint[0];

    data_y_pos = ScintPoint[1];
    data_y_pos_vuv = ScintPoint[1];
    data_y_pos_vis = ScintPoint[1];

    data_z_pos = ScintPoint[2];
    data_z_pos_vuv = ScintPoint[2];
    data_z_pos_vis = ScintPoint[2];

    data_event = event_number;
    data_event_vuv = event_number;
    data_event_vis = event_number;

	data_pmt = optical_channel;
	data_pmt_vuv = optical_channel;
	data_pmt_vis = optical_channel;

	// add entry for each photon
	for (int i = 0; i < num_VUV; i++) {
		data_time = times_vuv[i]; 
        data_time_vuv = times_vuv[i];
		data_tree_vuv->Fill();
        data_tree->Fill(); 
	}
	if (include_reflected) {
		for (int i = 0; i < num_VIS; i++) {
			data_time = times_vis[i]; 
            data_time_vis = times_vis[i];
			data_tree_vis->Fill();
	        data_tree->Fill(); 
		}
	}
}