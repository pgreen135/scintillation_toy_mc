#ifndef SOLID_ANGLE_H
#define SOLID_ANGLE_H

// class containing functions for calculating number of hits on each optical channel via the solid angle method
// calculates number of VUV hits, using gaisser-hillas corrections for Rayleigh scattering

// calculates number of visible hits via calculating number of hits on the cathode corrected with gaisser-hillas curves then the number
// of hits from the cathode for each optical channel using correction analogous to gaisser-hillas curves

#include <vector>
#include <string>

#include "TF1.h"
#include "TVector3.h"

class semi_analytic_hits {

private:

	// useful constants
	const double pi = 3.141592653589793;

	// detector type flag
	//const std::string flagDetector;

	// *************************************************************************************************
	//                    	OPTICAL DETECTOR SHAPE/SIZE
	// *************************************************************************************************	
	// Supercells: type = 1;
	double y_dimension_detector = 60;	// cm
	double z_dimension_detector = 60;	// cm
	// PMTs: type = 0; 
	double radius = 8*2.54/2.;	//	8" PMT diameter  // cm

	// structure definition for solid angle of rectangle function
	struct acc{
		// ax,ay,az = centre of rectangle; w = y dimension; h = z dimension
  		double ax, ay, az, w, h; 
	};

	bool _mathmore_loaded_ = false;


	// *************************************************************************************************
	//                    NUMBER OF VUV HITS PARAMETRIZATION
	// *************************************************************************************************
	// VUV hits Gaisser-Hillas Rayleigh scattering correction
	// Note: wires are now incorporated into the parameterisations, no longer need to add scaling factor
	
	// DUNE-VD Gaisser-Hillas
	// angle bins
	std::vector<double> angulo = {5, 15, 25, 35, 45, 55, 65, 75, 85};
	const double delta_angle = 10.;

	// Argon, RS = 99.9cm, flat PDs (Arapucas/Supercells)
	const double fGHVUVPars_flat_argon[4][9] = { {1.2343, 1.19482, 1.13939, 1.06642, 0.957333, 0.8411, 0.721859, 0.582155, 0.420655},
									{160.789, 163.041, 163.416, 176.419, 190.937, 205.741, 239.029, 244.506, 255.243}, 
									{90.6463, 92.4193, 94.1577, 139.613, 139.729, 188.447, 200, 200, 186.966}, 
									{-1000, -1000, -1000, -500, -500, -250, -100, -100, -100} }; 	
	const std::vector<double> slopes1_flat_argon = {-1.76996e-05, -1.53516e-05, -2.32847e-05, -2.15235e-05, -1.22799e-05, -2.12407e-05, -2.28983e-05, -1.17738e-05, -9.59168e-06};  
	const std::vector<double> slopes2_flat_argon = {-0.0342994, -0.0421957, -0.0253359, -0.0275521, -0.0354278, -0.0183091, -0.022059, -0.0228965, -0.0109844};  
	const std::vector<double> slopes3_flat_argon = {-0.0149828, -0.00686603, -0.0105492, -0.00897181, -0.00182785, -0.00985137, -3.68105e-07, 1.15155e-08, -0.00467523};  

    
    // Xenon, flat PDs (Arapucas/Supercells)
    const double fGHVUVPars_flat_xenon[4][9] = { {1.03634, 1.00335, 0.944009, 0.861764, 0.755889, 0.623379, 0.479978, 0.349492, 0.210756},
									{192.764, 179.4, 206.58, 227.04, 266.433, 345.513, 485.529, 634.413, 816.954}, 
									{70.0961, 84.8227, 82.3974, 158.585, 165.885, 183.76, 204.035, 200, 169.292},
									{-5000, -5000, -5000, -2500, -2500, -2500, -2500, -2500, -2500} };
	std::vector<double> slopes1_flat_xenon = {5.605e-05, 5.21965e-05, 4.23181e-05, 3.90189e-05, 3.34029e-05, 3.07638e-05, 1.89242e-05, 9.99064e-06, 1.52695e-06};  
	std::vector<double> slopes2_flat_xenon = {0.0901634, 0.12301, 0.0818468, 0.0712509, 0.0428058, 0.0448781, -0.014303, -0.0113832, 0.00367332};  
	std::vector<double> slopes3_flat_xenon = {-0.00829087, -0.0268017, -0.0123689, -0.0267884, -0.0189226, -0.0327574, -0.0103983, -8.18569e-07, 0.0310912};  							  

	// LAr absorption length in cm
	const double L_abs = 2000.;	// 20 m

	// *************************************************************************************************
	//                    NUMBER OF VIS HITS PARAMETRIZATION: N/A for VD
	// *************************************************************************************************
	// Note: wires are now incorporated into the parameterisations, no longer need to add scaling factor
	// Note: not applicable for the vertical drift
	
	// Detector properties: 
	
	// plane depth
	// Dune
	const double plane_depth = 363.38405;	// cm 
	
	// Cathode plane covered by foils
	// Dune
	// size
	const double y_dimension_foils = 1204.7255 + 5.466;	// cm		// 2 panels y height 602.36275 with 5.466cm gaps between
	const double z_dimension_foils = 1359.144 + 35.196;	// cm		// 6 panels of z width 226.524cm with 5.866cm gaps between them
	
	// centre coordinates
	// Dune
	const double x_foils = 363.38405; const double y_foils = 0; const double z_foils = 696.294;	// cm

	// DUNE-SP corrections
	// Argon, flat PDs (Arapucas/Supercells)
	const std::vector<double> vDistances_x_flat_argon = {0};			// cm	[16]
	const std::vector<double> vDistances_r_flat_argon = {0};		// cm	[9]

	const std::vector<std::vector<std::vector<double>>> fVISPars_flat_argon = {{ {0} }};

	// Xenon, flat PDs (Arapucas/Supercells)
	const std::vector<double> vDistances_x_flat_xenon = {0};			// cm	[16]
	const std::vector<double> vDistances_r_flat_xenon = {0};		// cm	[9]

	const std::vector<std::vector<std::vector<double>>> fVISPars_flat_xenon = {{ {0} }};
		
	

public:	
	// constructor 
	semi_analytic_hits();

	// destructor
	~semi_analytic_hits(){};

	// hits calculating functions
	int VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type, const int &optical_detector_orientation, const int &scintillation_type);
	int VisHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type, const int &scintillation_type);

	// gaisser-hillas function
	static Double_t GaisserHillas(double x, double *par);

	// solid angle of rectangular aperture calculation functions
	double omega(const double &a, const double &b, const double &d) const;
	double solid(const acc& out, const TVector3 &v) const;
	double solid_lateral(const acc& out, const TVector3 &v) const;

	// solid angle of circular aperture calculation functions
	double Disk_SolidAngle(double *x, double *p);
	double Disk_SolidAngle(double d, double h, double b);

	// solid angle of dome (PMTs)
	double Omega_Dome_Model(const double distance, const double theta) const;

	// linear interpolation function
	double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate );

};

#endif
