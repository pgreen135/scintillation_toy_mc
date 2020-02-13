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
	double y_dimension_detector = 9.3;	// cm
	double z_dimension_detector = 46.8;	// cm
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
	double fGHVUVPars[4][9];
	double fVUVBorderCorr[2];
	
	//  DUNE-SP Gaisser-Hillas
	// RS = 60cm							
	/*
	const double GH_SP[4][9] = { {1.37688, 1.35055, 1.29736, 1.21255, 1.0884, 0.985292, 0.828479, 0.686076, 0.531654}, 
										{105.721, 108.404, 114.289, 123.753, 142.937, 154.483, 179.057, 190.186, 188.029}, 
										{80, 80, 80, 80, 80, 80, 80, 80, 80}, 
										{-300, -300, -300, -300, -300, -300, -300, -300, -300} };
	*/
	// RS = 90cm
	
	const double GH_SP[4][9] = { {1.27459, 1.24136, 1.17434, 1.10557, 0.998743, 0.87692, 0.726141, 0.577982, 0.478801}, 
								  {111.684, 115.949, 139.661, 128.612, 144.167, 150.046, 186.602, 173.724, 125.866}, 
								  {45.7739, 46.0323, 41.8, 49.0639, 51.767, 59.6052, 59.5647, 75.981, 88.5809}, 
								  {-800, -800, -800, -800, -800, -800, -800, -800, -800} };
    
    // RS = 700cm
	/*							  
    const double GH_SP[4][9] = { {1.06944, 1.02935, 0.965427, 0.882084, 0.76177, 0.640781, 0.477896, 0.330519, 0.2039}, 
								  {42.5034, 74.6181, 101.596, 100.829, 142.166, 156.764, 248.771, 318.867, 308.229}, 
								  {135.609, 119.681, 100.001, 100.089, 100, 100.001, 100, 100.001, 100}, 
								  {-2500, -2500, -2500, -2500, -2500, -2500, -2500, -2500, -2500} };				
	*/
	
	// DUNE-SP border effects correction
	// RS = 60cm
	//const double BORDER_SP[2] = {9.58514e-05, 0.0858229};
	// RS =90cm
	const double BORDER_SP[2] = {6.31275e-05,  0.0482702};
	// RS = 700cm
	//const double BORDER_SP[2] = {-1.44353e-05,  -0.0157457};							  

	// bin size of the offset angle used in the parametrization in degrees
	const double delta_angulo = 10.;
	// LAr absorption length in cm
	const double L_abs = 2000.;

	// active volume corner
	const double fYactive_corner = 1200/2; 	// cm
	const double fZactive_corner = 1400/2;	// cm
	double fReference_to_corner;


	// *************************************************************************************************
	//                    NUMBER OF VIS HITS PARAMETRIZATION
	// *************************************************************************************************
	// properties of LAr for reflected light path calculation
	// refractive indices in LAr
	const double n_LAr_VUV = 2.632;     // effective index due to group vel.
	const double n_LAr_vis = 1.23;

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

	// Visible hits correction
	TF1* VIS_pol[9];
	const double delta_angle = 10.;

	// DUNE SP geometric correction [Arapuca/supercell, front window only]
	// RS = 60cm
	/*
	const double VIS_SP[6][9] = { {1.53686,1.46322,1.34083,1.14279,0.915688,0.648431,0.41744,0.41744,0.41744},
		{-0.0160284,-0.0145082,-0.0124793,-0.00975818,-0.0072178,-0.00464946,-0.00279529,-0.00279529,-0.00279529},
		{0.000112998,9.87981e-05,8.42566e-05,6.65798e-05,5.35068e-05,3.65524e-05,2.144e-05,2.144e-05,2.144e-05},
		{-4.75727e-07,-4.05306e-07,-3.49779e-07,-2.83337e-07,-2.52561e-07,-1.83123e-07,-1.04058e-07,-1.04058e-07,-1.04058e-07},
		{1.08826e-09,9.16466e-10,8.10421e-10,6.77421e-10,6.61672e-10,5.08153e-10,2.89974e-10,2.89974e-10,2.89974e-10},
		{-1.03393e-12,-8.73396e-13,-7.94945e-13,-6.86547e-13,-7.14112e-13,-5.73443e-13,-3.34853e-13,-3.34853e-13,-3.34853e-13} };
	*/
	// RS = 90cm
			
	const double VIS_SP[6][9] = { {1.48138, 1.44377, 1.34033, 1.14948, 0.930486, 0.671965, 0.439251, 0.439251, 0.439251},
		{-0.0129687, -0.0128759, -0.0122343, -0.00875131, -0.00635574, -0.00423979, -0.00268053, -0.00268053, -0.00268053},
		{8.49113e-05, 8.82379e-05, 9.62713e-05, 6.09814e-05, 4.93078e-05, 3.72279e-05, 2.50509e-05, 2.50509e-05, 2.50509e-05},
		{-3.61608e-07, -3.79819e-07, -4.76497e-07, -2.75142e-07, -2.52109e-07, -2.05943e-07, -1.4263e-07, -1.4263e-07, -1.4263e-07},
		{8.78898e-10, 9.03489e-10, 1.25367e-09, 6.96077e-10, 7.04725e-10, 6.04743e-10, 4.31868e-10, 4.31868e-10, 4.31868e-10},
		{-8.91232e-13, -8.88413e-13, -1.3072e-12, -7.2709e-13, -7.85031e-13, -6.92256e-13, -5.06101e-13, -5.06101e-13, -5.06101e-13} };
	
	// RS = 700cm
	/*	
	const double VIS_SP[6][9] = { {1.15911, 1.11916, 1.00025, 0.865688, 0.698012, 0.505619, 0.337874, 0.337874, 0.337874},
				{-0.00570387, -0.00588804, -0.00392452, -0.00267144, -0.00158738, -0.000358857, 0.000103533, 0.000103533, 0.000103533},
				{2.11059e-05, 3.35317e-05, 1.74076e-05, 1.60525e-05, 1.82463e-05, 8.42647e-06, 4.19648e-06, 4.19648e-06, 4.19648e-06},
				{-3.99733e-08, -1.4082e-07, -5.85585e-08, -8.2807e-08, -1.26157e-07, -6.25657e-08, -3.27901e-08, -3.27901e-08, -3.27901e-08},
				{2.10994e-11, 3.36018e-10, 1.22682e-10, 2.37932e-10, 3.99628e-10, 2.12365e-10, 1.24313e-10, 1.24313e-10, 1.24313e-10},
				{1.67885e-14, -3.27656e-13, -1.1701e-13, -2.68746e-13, -4.63135e-13, -2.65513e-13, -1.71392e-13, -1.71392e-13, -1.71392e-13} };
	*/	
	/*
	// ** 7m drift **
	// RS60cm
	const double VIS_RS60cm_SP[6][9] = { {1.40583,1.15359,1.25194,1.54649,1.70694,1.99575,1.99575,1.99575,1.99575},
	{-0.0073011,-0.00344252,-0.00411537,-0.00622707,-0.00588044,-0.00674468,-0.00674468,-0.00674468,-0.00674468},
	{3.21944e-05,1.0809e-05,1.46532e-05,2.42112e-05,2.14304e-05,2.87278e-05,2.87278e-05,2.87278e-05,2.87278e-05},
	{-9.17456e-08,-3.63564e-08,-4.78928e-08,-7.18351e-08,-6.77508e-08,-9.91708e-08,-9.91708e-08,-9.91708e-08,-9.91708e-08},
	{1.3461e-10,6.7504e-11,8.3892e-11,1.15271e-10,1.1713e-10,1.72574e-10,1.72574e-10,1.72574e-10,1.72574e-10},
	{-7.49764e-14,-4.42768e-14,-5.28968e-14,-6.92522e-14,-7.38904e-14,-1.074e-13,-1.074e-13,-1.074e-13,-1.074e-13} };	
	*/

	// DUNE SP border correction [Arapuca/supercell, front window only]
	// RS = 60 cm
	/*
	const std::vector<double> vDistances_x = {10, 40, 85, 130, 175, 220, 265, 310, 340, 355};			// cm	[10]
	const std::vector<double> vDistances_r = {122, 451, 647, 790};										// cm	[4]

	const std::vector<std::vector<std::vector<double>>> VIS_SP_Borders = { 
		{ {1.36687,1.22296,1.16675,1.14379,1.12301,1.10163,1.10099,1.09356,1.10187,1.12753},
		  {1.83254,1.61419,1.49286,1.43927,1.3912,1.35437,1.31543,1.27008,1.24884,1.25714},
		  {2.05931,1.85453,1.78544,1.79056,1.79707,1.8195,1.80777,1.7556,1.70645,1.70289},
		  {2.23968,2.15344,2.18786,2.36734,2.56407,2.75233,2.89681,2.90468,2.82028,2.79517}
		},
		{ {1.34339,1.20667,1.13506,1.11603,1.08991,1.07892,1.07119,1.06193,1.07188,1.09697},
		  {1.81072,1.59609,1.46559,1.4195,1.37173,1.32814,1.28887,1.2501,1.23295,1.24069},
		  {2.03827,1.82543,1.75276,1.74577,1.77848,1.78771,1.78426,1.73086,1.69444,1.68646},
		  {1.95378,1.84539,1.89026,2.02207,2.17497,2.33397,2.43329,2.42094,2.33427,2.31921}
		},
		{ {1.30212,1.1679,1.09268,1.06124,1.04809,1.02856,1.01896,1.00673,1.01583,1.03687},
		  {1.77623,1.53716,1.40967,1.36582,1.31276,1.28379,1.24822,1.21333,1.19201,1.19782},
		  {1.97939,1.75382,1.67763,1.70338,1.70174,1.71431,1.71404,1.67669,1.64612,1.64257},
		  {1.8617,1.72862,1.73262,1.86505,1.96976,2.09492,2.18476,2.16275,2.08759,2.06053}
		},
		{ {1.2217,1.09217,1.02034,1.0015,0.982541,0.969169,0.964904,0.960506,0.967135,0.985449},
		  {1.76632,1.55709,1.42782,1.38663,1.35201,1.32472,1.30128,1.2761,1.26203,1.26602},
		  {1.85359,1.6475,1.55838,1.56052,1.57037,1.57436,1.56395,1.53414,1.49935,1.49745},
		  {1.70973,1.60428,1.57908,1.68053,1.78143,1.87097,1.94755,1.94401,1.89018,1.874}
		},
		{ {1.19709,1.0421,0.975426,0.959783,0.943776,0.928057,0.925317,0.927927,0.939079,0.954253},
		  {1.80677,1.56675,1.45455,1.43069,1.40725,1.39217,1.38509,1.38628,1.39264,1.39857},
		  {1.78879,1.58522,1.50806,1.5289,1.54999,1.56731,1.58403,1.59104,1.58466,1.5856},
		  {1.51917,1.39067,1.3775,1.43522,1.52172,1.59911,1.65968,1.67423,1.63938,1.62868}
		},
		{ {1.31939,1.14616,1.06629,1.04966,1.03393,1.02082,1.0237,1.03905,1.05748,1.07774},
		  {1.44279,1.2267,1.11587,1.08317,1.0596,1.04128,1.03517,1.0428,1.05148,1.06037},
		  {1.49286,1.2971,1.20787,1.20624,1.20363,1.20469,1.20909,1.22173,1.2211,1.22693},
		  {1.35806,1.22853,1.19188,1.24884,1.30861,1.37408,1.43273,1.46065,1.45178,1.44975}
		},
		{ {1.67305,1.45578,1.35305,1.32855,1.32184,1.32226,1.33484,1.36954,1.40883,1.44865},
		  {1.539,1.30279,1.15584,1.11571,1.09565,1.08999,1.09669,1.11903,1.14344,1.16256},
		  {1.41221,1.20621,1.09857,1.08457,1.09274,1.10523,1.12587,1.15326,1.17331,1.18866},
		  {1.12174,0.991628,0.935007,0.961442,1.01713,1.07786,1.13943,1.18148,1.19394,1.2022}
		},
		{ {1.67305,1.45578,1.35305,1.32855,1.32184,1.32226,1.33484,1.36954,1.40883,1.44865},
		  {1.37135,1.13834,0.992489,0.951638,0.935147,0.933131,0.942681,0.970196,0.993699,1.01489},
		  {1.15366,0.96338,0.857661,0.841519,0.846182,0.86541,0.893007,0.92589,0.949376,0.968278},
		  {0.884385,0.759848,0.707666,0.724935,0.767402,0.824071,0.883741,0.935287,0.958996,0.97325}
		},
		{ {1.67305,1.45578,1.35305,1.32855,1.32184,1.32226,1.33484,1.36954,1.40883,1.44865},
		  {1.37135,1.13834,0.992489,0.951638,0.935147,0.933131,0.942681,0.970196,0.993699,1.01489},
		  {1.15366,0.96338,0.857661,0.841519,0.846182,0.86541,0.893007,0.92589,0.949376,0.968278},
		  {0.884385,0.759848,0.707666,0.724935,0.767402,0.824071,0.883741,0.935287,0.958996,0.97325}
		}
	};
	*/
	
	// RS 90cm
	
	const std::vector<double> vDistances_x = {10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350};	// cm	[15]
	const std::vector<double> vDistances_r = {123, 453, 646, 793};			// cm	[4]

	const std::vector<std::vector<std::vector<double>>> VIS_SP_Borders = { 
		{ {1.38887, 1.28656, 1.22529, 1.19441, 1.17094, 1.16426, 1.14155, 1.12663, 1.11825, 1.1028, 1.0997, 1.09631, 1.09442, 1.09982, 1.12772},
		  {1.86426, 1.71695, 1.60912, 1.5239, 1.49921, 1.45715, 1.44296, 1.41008, 1.3696, 1.35278, 1.3124, 1.2966, 1.27699, 1.2675, 1.2712},
		  {2.10302, 2.00165, 1.92391, 1.89817, 1.89424, 1.91369, 1.91249, 1.94263, 1.94218, 1.92922, 1.93573, 1.93887, 1.92529, 1.90906, 1.88003},
		  {2.54832, 2.46558, 2.43422, 2.41528, 2.52735, 2.62071, 2.69887, 2.76163, 2.82774, 2.89572, 2.95643, 2.96838, 2.92766, 2.88496, 2.8223}
		},
		{ {1.36558, 1.28038, 1.20765, 1.17506, 1.14436, 1.13632, 1.11304, 1.10529, 1.07668, 1.0756, 1.06931, 1.07555, 1.06838, 1.07234, 1.08735},
		  {1.83774, 1.7037, 1.54094, 1.50718, 1.45678, 1.42634, 1.39743, 1.36315, 1.32885, 1.32497, 1.30464, 1.29334, 1.26518, 1.24804, 1.23667},
		  {2.02689, 1.92889, 1.81812, 1.7748, 1.77256, 1.77907, 1.7718, 1.77088, 1.7585, 1.77364, 1.75921, 1.75915, 1.74151, 1.70829, 1.67404},
		  {2.19339, 2.16153, 2.15268, 2.12713, 2.14764, 2.20089, 2.26327, 2.3332, 2.36635, 2.44178, 2.47865, 2.48745, 2.48303, 2.40155, 2.30349}
		},
		{ {1.30385, 1.2129, 1.12647, 1.11346, 1.09402, 1.08805, 1.07607, 1.05912, 1.02795, 1.02837, 1.02569, 1.02581, 1.02908, 1.02919, 1.02763},
		  {1.75594, 1.59184, 1.48457, 1.44582, 1.40791, 1.38835, 1.36028, 1.32781, 1.29448, 1.28538, 1.26454, 1.2573, 1.23627, 1.2162, 1.19196},
		  {1.94354, 1.79952, 1.70626, 1.66641, 1.6887, 1.68805, 1.68825, 1.67636, 1.67394, 1.67826, 1.67358, 1.6833, 1.66695, 1.62526, 1.57446},
		  {2.09017, 1.94979, 1.92513, 1.90346, 1.9645, 2.04143, 2.07791, 2.12547, 2.15142, 2.21473, 2.23871, 2.27191, 2.23392, 2.17906, 2.06693}
		},
		{ {1.21911, 1.15514, 1.08586, 1.05224, 1.03327, 1.01892, 1.0038, 0.989274, 0.980393, 0.971908, 0.971483, 0.969728, 0.967443, 0.964384, 0.970291},
		  {1.74303, 1.63335, 1.52126, 1.46386, 1.42366, 1.40107, 1.38285, 1.35969, 1.33874, 1.33097, 1.30987, 1.29225, 1.28457, 1.26475, 1.25129},
		  {1.82069, 1.7568, 1.66142, 1.62711, 1.62267, 1.62074, 1.6323, 1.63966, 1.64853, 1.65429, 1.66901, 1.6648, 1.65679, 1.62547, 1.5876},
		  {1.885, 1.82586, 1.76656, 1.77214, 1.79626, 1.84132, 1.87874, 1.91991, 1.95655, 1.96812, 1.99918, 1.9967, 1.97199, 1.91895, 1.84386}
		},
		{ {1.20502, 1.12323, 1.04705, 1.00297, 0.988905, 0.977234, 0.967661, 0.952328, 0.941662, 0.933208, 0.927596, 0.932176, 0.929665, 0.931696, 0.938101},
		  {1.77596, 1.66357, 1.53535, 1.4664, 1.44399, 1.4304, 1.42068, 1.40318, 1.38439, 1.36726, 1.37582, 1.37428, 1.3745, 1.37521, 1.36774},
		  {1.68525, 1.59554, 1.49654, 1.46611, 1.46066, 1.46367, 1.46301, 1.47848, 1.47297, 1.48265, 1.48248, 1.49203, 1.49091, 1.47455, 1.43938},
		  {1.67028, 1.56903, 1.52651, 1.53562, 1.53816, 1.56908, 1.61138, 1.64199, 1.65008, 1.66909, 1.70485, 1.69958, 1.69162, 1.64988, 1.58941}
		},
		{ {1.30455, 1.21736, 1.12749, 1.09126, 1.07541, 1.06265, 1.05353, 1.04232, 1.03417, 1.02454, 1.02441, 1.02724, 1.03628, 1.04361, 1.05003},
		  {1.39963, 1.30173, 1.19587, 1.14253, 1.11904, 1.08834, 1.07596, 1.05846, 1.05212, 1.03273, 1.03211, 1.02958, 1.03722, 1.03737, 1.03211},
		  {1.42469, 1.3497, 1.25591, 1.21902, 1.21221, 1.20994, 1.21021, 1.20938, 1.20737, 1.20433, 1.21053, 1.21232, 1.21929, 1.21327, 1.19146},
		  {1.49179, 1.41844, 1.33146, 1.32161, 1.32748, 1.36183, 1.37731, 1.40147, 1.41801, 1.42756, 1.44585, 1.46643, 1.45562, 1.43464, 1.39811}
		},
		{ {1.6401, 1.52713, 1.41132, 1.36182, 1.3366, 1.33563, 1.33416, 1.32711, 1.32154, 1.31632, 1.3248, 1.34096, 1.35568, 1.37485, 1.39351},
		  {1.49175, 1.36917, 1.2463, 1.17668, 1.14579, 1.12082, 1.11192, 1.09797, 1.08352, 1.07673, 1.08222, 1.08712, 1.10306, 1.11773, 1.11685},
		  {1.32827, 1.24325, 1.13731, 1.08398, 1.07668, 1.06424, 1.0718, 1.0704, 1.06985, 1.07835, 1.08543, 1.09874, 1.11195, 1.11573, 1.1029},
		  {1.21849, 1.15967, 1.07397, 1.0455, 1.04265, 1.04644, 1.07214, 1.08678, 1.10662, 1.12403, 1.14751, 1.15754, 1.17154, 1.16163, 1.14496}
		},
		{ {1.6401, 1.52713, 1.41132, 1.36182, 1.3366, 1.33563, 1.33416, 1.32711, 1.32154, 1.31632, 1.3248, 1.34096, 1.35568, 1.37485, 1.39351},
		  {1.33359, 1.21674, 1.08101, 1.01555, 0.982123, 0.961531, 0.961298, 0.936583, 0.925483, 0.923806, 0.921805, 0.93999, 0.951694, 0.95851, 0.977982},
		  {1.12958, 1.05089, 0.944507, 0.891667, 0.878727, 0.868075, 0.865406, 0.872593, 0.875672, 0.884941, 0.897938, 0.91557, 0.935379, 0.942483, 0.948206},
		  {0.96415, 0.907609, 0.82139, 0.791617, 0.788006, 0.788812, 0.804325, 0.819716, 0.837583, 0.859831, 0.878327, 0.897308, 0.921485, 0.926272, 0.919751}
		},
		{ {1.6401, 1.52713, 1.41132, 1.36182, 1.3366, 1.33563, 1.33416, 1.32711, 1.32154, 1.31632, 1.3248, 1.34096, 1.35568, 1.37485, 1.39351},
		  {1.33359, 1.21674, 1.08101, 1.01555, 0.982123, 0.961531, 0.961298, 0.936583, 0.925483, 0.923806, 0.921805, 0.93999, 0.951694, 0.95851, 0.977982},
		  {1.12958, 1.05089, 0.944507, 0.891667, 0.878727, 0.868075, 0.865406, 0.872593, 0.875672, 0.884941, 0.897938, 0.91557, 0.935379, 0.942483, 0.948206},
		  {0.96415, 0.907609, 0.82139, 0.791617, 0.788006, 0.788812, 0.804325, 0.819716, 0.837583, 0.859831, 0.878327, 0.897308, 0.921485, 0.926272, 0.919751}
		}
	};
	

	// RS = 700 cm
	/*
	const std::vector<double> vDistances_x = {10, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350};	// cm	[15]
	const std::vector<double> vDistances_r = {123, 453, 646, 793};			// cm	[4]

	const std::vector<std::vector<std::vector<double>>> VIS_SP_Borders = { 
		{ {1.34257, 1.33392, 1.30955, 1.26805, 1.24494, 1.21715, 1.19967, 1.17366, 1.16278, 1.14946, 1.13491, 1.12195, 1.11507, 1.11107, 1.13155},
		  {1.81434, 1.80089, 1.74854, 1.68386, 1.62008, 1.56887, 1.5142, 1.49742, 1.4546, 1.4148, 1.38732, 1.34481, 1.30965, 1.27681, 1.26784},
		  {2.10302, 2.00165, 1.92391, 1.89817, 1.89424, 1.91369, 1.91249, 1.94263, 1.94218, 1.92922, 1.93573, 1.93887, 1.92529, 1.90906, 1.88003},
		  {2.54832, 2.46558, 2.43422, 2.41528, 2.52735, 2.62071, 2.69887, 2.76163, 2.82774, 2.89572, 2.95643, 2.96838, 2.92766, 2.88496, 2.8223}
		},
		{ {1.32771, 1.31711, 1.27665, 1.23299, 1.21882, 1.19228, 1.17606, 1.14836, 1.13515, 1.11743, 1.10733, 1.09598, 1.09313, 1.08534, 1.09584},
		  {1.80206, 1.77589, 1.70169, 1.62308, 1.59778, 1.56933, 1.51024, 1.4693, 1.42141, 1.38222, 1.34837, 1.3232, 1.30592, 1.26406, 1.23888},
		  {2.02689, 1.92889, 1.81812, 1.7748, 1.77256, 1.77907, 1.7718, 1.77088, 1.7585, 1.77364, 1.75921, 1.75915, 1.74151, 1.70829, 1.67404},
		  {2.19339, 2.16153, 2.15268, 2.12713, 2.14764, 2.20089, 2.26327, 2.3332, 2.36635, 2.44178, 2.47865, 2.48745, 2.48303, 2.40155, 2.30349}
		},
		{ {1.28848, 1.2626, 1.21977, 1.18287, 1.15554, 1.13857, 1.11822, 1.10205, 1.09142, 1.07503, 1.06215, 1.04743, 1.04262, 1.04084, 1.05145},
		  {1.72437, 1.68789, 1.62692, 1.57661, 1.51298, 1.46886, 1.4456, 1.40796, 1.38129, 1.3606, 1.31937, 1.29062, 1.2682, 1.2419, 1.21801},
		  {1.94354, 1.79952, 1.70626, 1.66641, 1.6887, 1.68805, 1.68825, 1.67636, 1.67394, 1.67826, 1.67358, 1.6833, 1.66695, 1.62526, 1.57446},
		  {0.09017, 1.94979, 1.92513, 1.90346, 1.9645, 2.04143, 2.07791, 2.12547, 2.15142, 2.21473, 2.23871, 2.27191, 2.23392, 2.17906, 2.06693}
		},
		{ {1.22651, 1.20078, 1.15347, 1.12557, 1.09833, 1.08271, 1.0633, 1.0467, 1.03215, 1.01808, 1.00657, 0.99631, 0.990112, 0.988301, 0.994004},
		  {1.73689, 1.69324, 1.63096, 1.58067, 1.54511, 1.51551, 1.48217, 1.45228, 1.41523, 1.38207, 1.36127, 1.33143, 1.31365, 1.29692, 1.27638},
		  {1.82069, 1.7568, 1.66142, 1.62711, 1.62267, 1.62074, 1.6323, 1.63966, 1.64853, 1.65429, 1.66901, 1.6648, 1.65679, 1.62547, 1.5876},
		  {1.885, 1.82586, 1.76656, 1.77214, 1.79626, 1.84132, 1.87874, 1.91991, 1.95655, 1.96812, 1.99918, 1.9967, 1.97199, 1.91895, 1.84386}
		},
		{ {1.20408, 1.1635, 1.10885, 1.0857, 1.06482, 1.04906, 1.0338, 1.01726, 1.00101, 0.984229, 0.9713, 0.969661, 0.964979, 0.968915, 0.975611},
		  {1.72266, 1.68792, 1.62566, 1.58473, 1.5516, 1.53852, 1.51755, 1.488, 1.4631, 1.4446, 1.42138, 1.41905, 1.42276, 1.40947, 1.40578},
		  {1.68525, 1.59554, 1.49654, 1.46611, 1.46066, 1.46367, 1.46301, 1.47848, 1.47297, 1.48265, 1.48248, 1.49203, 1.49091, 1.47455, 1.43938},
		  {1.67028, 1.56903, 1.52651, 1.53562, 1.53816, 1.56908, 1.61138, 1.64199, 1.65008, 1.66909, 1.70485, 1.69958, 1.69162, 1.64988, 1.58941}
		},
		{ {1.29399, 1.25973, 1.2065, 1.17912, 1.15625, 1.13997, 1.1308, 1.12112, 1.10293, 1.09608, 1.09075, 1.08472, 1.08388, 1.08806, 1.10483},
		  {1.36599, 1.34403, 1.29717, 1.24754, 1.21609, 1.19433, 1.17798, 1.15062, 1.12842, 1.11128, 1.09633, 1.08446, 1.08041, 1.0734, 1.07428},
		  {1.42469, 1.3497, 1.25591, 1.21902, 1.21221, 1.20994, 1.21021, 1.20938, 1.20737, 1.20433, 1.21053, 1.21232, 1.21929, 1.21327, 1.19146},
		  {1.49179, 1.41844, 1.33146, 1.32161, 1.32748, 1.36183, 1.37731, 1.40147, 1.41801, 1.42756, 1.44585, 1.46643, 1.45562, 1.43464, 1.39811}
		},
		{ {1.56023, 1.52742, 1.46551, 1.43774, 1.41781, 1.41213, 1.40432, 1.40137, 1.40054, 1.39587, 1.39967, 1.40396, 1.41416, 1.43828, 1.47772},
		  {1.45098, 1.42689, 1.36724, 1.32312, 1.28266, 1.25574, 1.23895, 1.22242, 1.19771, 1.18354, 1.17823, 1.17356, 1.16905, 1.18034, 1.19853},
		  {1.32827, 1.24325, 1.13731, 1.08398, 1.07668, 1.06424, 1.0718, 1.0704, 1.06985, 1.07835, 1.08543, 1.09874, 1.11195, 1.11573, 1.1029},
		  {1.21849, 1.15967, 1.07397, 1.0455, 1.04265, 1.04644, 1.07214, 1.08678, 1.10662, 1.12403, 1.14751, 1.15754, 1.17154, 1.16163, 1.14496}
		},
		{ {1.56023, 1.52742, 1.46551, 1.43774, 1.41781, 1.41213, 1.40432, 1.40137, 1.40054, 1.39587, 1.39967, 1.40396, 1.41416, 1.43828, 1.47772},
		  {1.35078, 1.33182, 1.24785, 1.19712, 1.16909, 1.13515, 1.11397, 1.08589, 1.07972, 1.06323, 1.05141, 1.0467, 1.03599, 1.0631, 1.07381},
		  {1.12958, 1.05089, 0.944507, 0.891667, 0.878727, 0.868075, 0.865406, 0.872593, 0.875672, 0.884941, 0.897938, 0.91557, 0.935379, 0.942483, 0.948206},
		  {0.96415, 0.907609, 0.82139, 0.791617, 0.788006, 0.788812, 0.804325, 0.819716, 0.837583, 0.859831, 0.878327, 0.897308, 0.921485, 0.926272, 0.919751}
		},
		{ {1.56023, 1.52742, 1.46551, 1.43774, 1.41781, 1.41213, 1.40432, 1.40137, 1.40054, 1.39587, 1.39967, 1.40396, 1.41416, 1.43828, 1.47772},
		  {1.35078, 1.33182, 1.24785, 1.19712, 1.16909, 1.13515, 1.11397, 1.08589, 1.07972, 1.06323, 1.05141, 1.0467, 1.03599, 1.0631, 1.07381},
		  {0.12958, 1.05089, 0.944507, 0.891667, 0.878727, 0.868075, 0.865406, 0.872593, 0.875672, 0.884941, 0.897938, 0.91557, 0.935379, 0.942483, 0.948206},
		  {0.96415, 0.907609, 0.82139, 0.791617, 0.788006, 0.788812, 0.804325, 0.819716, 0.837583, 0.859831, 0.878327, 0.897308, 0.921485, 0.926272, 0.919751}
		} 
	};
	*/	
	

public:	
	// constructor 
	semi_analytic_hits();

	// destructor
	~semi_analytic_hits(){};

	// hits calculating functions
	int VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type);
	int VisHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type);

	// gaisser-hillas function
	static Double_t GaisserHillas(double x, double *par);

	// solid angle of rectangular aperture calculation functions
	double omega(const double &a, const double &b, const double &d) const;
	double solid(const acc& out, const TVector3 &v) const;

	// solid angle of circular aperture calculation functions
	double Disk_SolidAngle(double *x, double *p);
	double Disk_SolidAngle(double d, double h, double b);

	// linear interpolation function
	double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate );

};

#endif