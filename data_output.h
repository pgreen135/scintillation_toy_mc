#ifndef DATA_OUTPUT_H
#define DATA_OUTPUT_H

// class handling writing of event information and detected photon information to output root file

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

class data_output {

private:
	// output file
	TFile *output_file;
	
	// branches:
	TTree *event_tree;
	TTree *data_tree;
	TTree *data_tree_vuv;
	TTree *data_tree_vis;

	// tree entries:
	// event
	int event_no;
	double event_E; 
	double event_time;
	int event_particle_type;
	double event_x_pos, event_y_pos, event_z_pos;
	// data
	int data_pmt, data_pmt_vuv, data_pmt_vis;
	int data_event, data_event_vuv, data_event_vis;
	double data_time, data_time_vuv, data_time_vis;
	double data_x_pos, data_x_pos_vuv, data_x_pos_vis;
	double data_y_pos, data_y_pos_vuv, data_y_pos_vis;
	double data_z_pos, data_z_pos_vuv, data_z_pos_vis;

	const bool include_reflected;

public:

	// constructor
	data_output(const char* output_file_name, const bool include_timings, const bool include_reflected);

	// destructor 
	~data_output();

	// add event tree entry
	void add_event(const int &event_number, const double &event_energy, const double &time, const double &particle_type, const std::vector<double> &event_position);

	// add data tree entry
	// without times
	void add_data(const int &event_number, const int &optical_channel, const int &num_VUV, const int &num_VIS, const TVector3 &ScintPoint);
	// with times
	void add_data(const int &event_number, const int &optical_channel, const int &num_VUV, const int &num_VIS, const TVector3 &ScintPoint, const std::vector<double> &times_vuv, const std::vector<double> &times_vis);

	// write file
	void write_output_file() { output_file->Write(); }

};
















#endif