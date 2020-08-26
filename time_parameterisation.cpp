// implementation of time parameterisation class

#include "time_parameterisation.h"
#include "time_parameterisation_services.h"

#include <string>
#include <cmath>

#include "TMath.h"
#include "TRandom.h"

using namespace std;

// constructor
time_parameterisation::time_parameterisation(const double &size): step_size{size} {	
	
    // create vector of empty TF1s that will be replaces with the parameterisations that are generated as they are required
    // default TF1() constructor gives function with 0 dimensions, can then check numDim to qucikly see if a parameterisation has been generated  
    int num_params = (d_max - 25) / step_size;  // for d < 25cm, no parameterisaton, a delta function is used instead	
    vector<TF1> VUV_timing_temp(num_params,TF1());
	// insert empty TF1 vectors, for each angular bin
    for (int i = 0; i < timing_vuv_angular_bins; i++) {
        VUV_timing.push_back(VUV_timing_temp);    
    }   
    
    // initialise vectors to contain range parameterisations sampled to in each case
    // when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling is regenerated, this is the slow part!
    vector<double> VUV_empty(num_params, 0);
    for (int i = 0; i < timing_vuv_angular_bins; i++) {
        VUV_max.push_back(VUV_empty);
        VUV_min.push_back(VUV_empty);
    }
}

// parameterisation generation function
void time_parameterisation::generateparam(const int &index, const int &angle_bin) {
   
    // get distance 
    double distance_in_cm = (index * step_size) + 25;
    
    // time range
    const double signal_t_range = 5000.;
   
    // parameterisation TF1    
    TF1 fVUVTiming;
      
    // For very short distances the time correction is just a shift
    double t_direct_mean = distance_in_cm/vuv_vgroup_mean;
    double t_direct_min = distance_in_cm/vuv_vgroup_max;
      
    // Defining the model function(s) describing the photon transportation timing vs distance 
    // Getting the landau parameters from the time parametrization
    double* pars_landau = interpolate(vDistances_landau, vMpv_landau[angle_bin], vWidth_landau[angle_bin], vNorm_landau[angle_bin], distance_in_cm, true);
    // Deciding which time model to use (depends on the distance)
    // defining useful times for the VUV arrival time shapes
    if(distance_in_cm >= inflexion_point_distance) {
        double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
        // Set model: Landau 
        fVUVTiming =  TF1("fVUVTiming",model_far,0,signal_t_range,4);
        fVUVTiming.SetParameters(pars_far);
    }
    else {
        // Set model: Landau + Exponential 
        fVUVTiming =  TF1("fVUVTiming",model_close,0,signal_t_range,7); 
        // Exponential parameters
        double pars_expo[2];   
        // Getting the exponential parameters from the time parametrization
        pars_expo[1] = interpolate(vDistances_exp, vSlope_exp[angle_bin], distance_in_cm, true);
        pars_expo[0] = interpolate(vDistances_exp, vNorm_exp[angle_bin], distance_in_cm, true);
        pars_expo[0] *= pars_landau[2];
        pars_expo[0] = log(pars_expo[0]);
        // this is to find the intersection point between the two functions:
        TF1 fint = TF1("fint",finter_d,pars_landau[0],4*t_direct_mean,5);
        double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
        fint.SetParameters(parsInt);
        double t_int = fint.GetMinimumX();
        double minVal = fint.Eval(t_int);
        // the functions must intersect - output warning if they don't
        if(minVal>0.015) {
            std::cout<<"WARNING: Parametrization of VUV light discontinuous for distance = " << distance_in_cm << std::endl;
        }
        //delete fint;   
        double parsfinal[7] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1], t_direct_min};
        fVUVTiming.SetParameters(parsfinal);    
    }

    // calculate max and min distance relevant to sample parameterisation 
    // max
    const int nq_max=1;
    double xq_max[nq_max];
    double yq_max[nq_max];    
    xq_max[0] = 0.995;   // include 99.5% of tail
    fVUVTiming.GetQuantiles(nq_max,yq_max,xq_max);
    double max = yq_max[0];
    // min
    double min = t_direct_min;

    // set the number of points used to sample parameterisation
    // for shorter distances, peak is sharper so more sensitive sampling required - values could be optimised, but since these are only generate once difference is not significant
    int f_sampling;
    if (distance_in_cm < 50) { f_sampling = 10000; }
    else if (distance_in_cm < 100){ f_sampling = 5000; }
    else { f_sampling = 1000; }
    fVUVTiming.SetNpx(f_sampling);    

    // generate the sampling
    // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
    // all subsequent calls check if it has been generated previously and are ~100+ times quicker
    double arrival_time = fVUVTiming.GetRandom(min,max);
    // add timing to the vector of timings and range to vectors of ranges
    VUV_timing[angle_bin][index] = fVUVTiming;
    VUV_max[angle_bin][index] = max;
    VUV_min[angle_bin][index] = min;
    
}

// VUV arrival times calculation function
vector<double> time_parameterisation::getVUVTime(const double &distance, const int &angle_bin, const int &number_photons) {
    
    // pre-allocate memory
    std::vector<double> arrival_time_distrb;
    arrival_time_distrb.clear();
    arrival_time_distrb.reserve(number_photons);

    // distance < 25cm
    if (distance < 25) {
        // times are fixed shift i.e. direct path only
        double t_prop_correction = distance/vuv_vgroup_mean;
        for (int i = 0; i < number_photons; i++){
            arrival_time_distrb.push_back(t_prop_correction);
        }
    }
    // distance >= 25cm
    else {
        // determine nearest parameterisation in discretisation
        int index = std::round((distance - 25) / step_size);
        // check whether required parameterisation has been generated, generating if not
        if (VUV_timing[angle_bin][index].GetNdim() == 0) {
            generateparam(index, angle_bin);
        }
        // randomly sample parameterisation for each photon
        for (int i = 0; i < number_photons; i++){
            arrival_time_distrb.push_back(VUV_timing[angle_bin][index].GetRandom(VUV_min[angle_bin][index],VUV_max[angle_bin][index]));
        }  
    }
    return arrival_time_distrb;
}

// vis arrival times calculation function
vector<double> time_parameterisation::getVisTime(const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &number_photons) {
    
    // *************************************************************************************************
    // Calculation of earliest arrival times and corresponding unsmeared distribution
    // *************************************************************************************************

    // bounce_point is point on wall where TPB is activated most intensely by the scintillation
    TVector3 bounce_point(cathode_plane_depth,ScintPoint[1],ScintPoint[2]);
   
    // calculate distance travelled by VUV light and by vis light
    double VUVdist = (bounce_point-ScintPoint).Mag();
    double Visdist = (OpDetPoint-bounce_point).Mag();

    // calculate times taken by each part
    int angle_bin_vuv = 0;   // on-axis by definition
    vector<double> VUVTimes  = getVUVTime(VUVdist, angle_bin_vuv, number_photons);
    vector<double> ReflTimes(number_photons, Visdist/vis_vmean);

    // sum parts to get total transport times times                    
    vector<double> transport_time_vis(number_photons, 0.0);
    for (int i=0; i<number_photons; i++) {
        transport_time_vis[i] = VUVTimes[i] + ReflTimes[i];
    }
    

    // *************************************************************************************************
    // Smearing of arrival time distribution
    // *************************************************************************************************
  
    // calculate fastest time possible
    // vis part
    double vis_time = Visdist/vis_vmean;
    // vuv part
    double vuv_time;
    if (VUVdist < 25){
        vuv_time = VUVdist/vuv_vgroup_max;
    }
    else {
        // find index of required parameterisation
        int index = std::round((VUVdist - 25) / step_size);
        // find shortest time
        vuv_time = getVUVmin(index, angle_bin_vuv);
    }
    // sum
    double fastest_time = vis_time + vuv_time;

    // calculate angle theta between reflection point and optical detector
    double cosine_theta = std::abs(OpDetPoint[0] - bounce_point[0]) / Visdist;
    double theta = acos(cosine_theta)*180./pi;
    
    // determine smearing parameters using interpolation of generated points: 
    //  1). tau = exponential smearing factor, varies with distance and angle
    //  2). cutoff = largest smeared time allowed, preventing excessively large times caused by exponential
    // distance to cathode
    double distance_cathode_plane = std::abs(cathode_plane_depth - ScintPoint[0]);
    // angular bin
    int theta_bin = theta / timing_vis_angular_bin_size;
    // radial distance from centre of TPC (y,z plane)
    double R = sqrt( pow(ScintPoint[1],2) + pow(ScintPoint[2] - z_foils, 2) );

    // cut-off and tau
    // cut-off
    // interpolate in x for each r bin
    std::vector<double> interp_vals(cut_off_bins[theta_bin].size(), 0.0);
    for (int i = 0; i < cut_off_bins[theta_bin].size(); i++){
        interp_vals[i] = interpolate(refl_vdistances, cut_off_bins[theta_bin][i], distance_cathode_plane, true);
    }
    // interpolate in r
    double cutoff = interpolate(refl_vradial, interp_vals, R, true);
    
    // tau
    // interpolate in x for each r bin
    std::vector<double> interp_vals_tau(tau_bins[theta_bin].size(), 0.0);
    for (int i = 0; i < tau_bins[theta_bin].size(); i++){
        interp_vals_tau[i] = interpolate(refl_vdistances, tau_bins[theta_bin][i], distance_cathode_plane, true);
    }
    // interpolate in r
    double tau = interpolate(refl_vradial, interp_vals_tau, R, true);

    // apply smearing:
    for (int i=0; i < number_photons; i++){
        double arrival_time = transport_time_vis[i];
        double arrival_time_smeared;
        // if time is already greater than cutoff or minimum smeared time would be greater than cutoff, do not apply smearing
        if (arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(1.0))-1) >= cutoff) {
            arrival_time_smeared = arrival_time;
        }
        // otherwise smear
        else {
            int counter = 0;
            // loop until time generated is within cutoff limit
            // most are within single attempt, very few take more than two -- could be made more efficient, some way of avoiding sluggish do-while loop
            do {
                // don't attempt smearings too many times for cases near cutoff (very few cases, not smearing these makes negigible difference)
                if (counter >= 10){
                    arrival_time_smeared = arrival_time; // don't smear
                    break;
                }
                else {
                    // generate random number in appropriate range            
                    double x = gRandom->Uniform(0.5,1.0);
                    // apply the exponential smearing
                    arrival_time_smeared = arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(x))-1);
                }
                // increment counter
                counter++;
            }   while (arrival_time_smeared > cutoff);
        }
        transport_time_vis[i] = arrival_time_smeared;
    }  
    
    return transport_time_vis;
}

double time_parameterisation::getVUVmin(const int &index, const int &angle_bin){
    
    if (VUV_min[angle_bin][index] == 0) {
        generateparam(index, angle_bin);
    }   
    double min = VUV_min[angle_bin][index];
    return min;
}
