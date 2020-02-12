#include "semi_analytic_hits.h"

// implementation of semi-analytic model for number of incident photons

#include <iostream>
#include <cmath>

#include "TRandom.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFormula.h"
#include "Math/SpecFuncMathMore.h"

using namespace std;

// constructor
semi_analytic_hits::semi_analytic_hits() {

  // load mathmore library
  gSystem->Load("libMathMore.so");
  if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
  _mathmore_loaded_ = true;

  std::cout<<"Light simulation for DUNE Single Phase detector"<<std::endl;

  // initialise gaisser hillas functions for VUV Rayleigh scattering correction
  // DUNE-SP
  fReference_to_corner = sqrt(pow(fYactive_corner, 2) + pow(fZactive_corner, 2));
  for (int bin = 0; bin < 9; bin ++) {
    for (int j = 0; j < 4; j++) {
      fGHVUVPars[j][bin] = GH_RS60cm_SP[j][bin];
    }
  }
  for (int i = 0; i < 2; i++) {
    fVUVBorderCorr[i] = BORDER_RS60cm_SP[i];
  }
  
  // initialise pol5 functions for visible hits correction
  double pars_ini_vis[6] = {0,0,0,0,0,0};
  for (int bin = 0; bin < 9; bin++) {
    VIS_pol[bin] = new TF1 ("pol", "pol5", 0, 2000);
    for (int j = 0; j < 6; j++){
      pars_ini_vis[j] = VIS_RS60cm_SP[j][bin];
    }
    VIS_pol[bin]->SetParameters(pars_ini_vis);
  }

  std::cout << std::endl;

}

// VUV hits calculation
int semi_analytic_hits::VUVHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type) {
  
  // distance and angle between ScintPoint and OpDetPoint
  double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
  double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance;
  double theta = acos(cosine)*180./pi;

  // calculate solid angle:
  double solid_angle = 0;
  // rectangular aperture
  if (optical_detector_type == 1) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint; 
    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
    detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window

    // get scintillation point coordinates relative to arapuca window centre
    TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;  

    // calculate solid angle
    solid_angle = solid(detPoint, ScintPoint_rel);
  }
  // disk aperture
  else if (optical_detector_type == 0) {
    // offset in z-y plane
    double d = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
    // drift distance (in x)
    double h =  sqrt(pow(ScintPoint[0] - OpDetPoint[0],2));
    // Solid angle of a disk
    solid_angle = Disk_SolidAngle(d, h, radius);
  }
  else {
    std::cout << "Error: Invalid optical detector type." << endl;
    exit(1);
  }  

  // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
  double hits_geo = exp(-1.*distance/L_abs) * (solid_angle / (4*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = (theta/delta_angulo);
  
  // account for border effects in Gaisser-Hillas function
  double z_to_corner = abs(ScintPoint[2] - fZactive_corner) - fZactive_corner;
  double y_to_corner = abs(ScintPoint[1]) - fYactive_corner;
  double distance_to_corner = sqrt(y_to_corner*y_to_corner + z_to_corner*z_to_corner);  // in the ph-cathode plane
  
  double pars_ini_[4] = {fGHVUVPars[0][j] + fVUVBorderCorr[0] * (distance_to_corner - fReference_to_corner),
                         fGHVUVPars[1][j] + fVUVBorderCorr[1] * (distance_to_corner - fReference_to_corner),
                         fGHVUVPars[2][j],
                         fGHVUVPars[3][j]};
  
  double GH_correction = GaisserHillas(distance, pars_ini_);
  
  // apply correction
  double hits_rec = gRandom->Poisson( GH_correction*hits_geo/cosine );

  // round to integer value, cannot have non-integer number of hits
  int hits_vuv = std::round(hits_rec);

  return hits_vuv;
}

// Visible hits calculation
int semi_analytic_hits::VisHits(const int &Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint, const int &optical_detector_type) {
 
  // 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:

  // set cathode plane struct for solid angle function
  acc cathode_plane; 
  cathode_plane.ax = x_foils; cathode_plane.ay = y_foils; cathode_plane.az = z_foils;   		// centre coordinates of cathode plane
  cathode_plane.w = y_dimension_foils; cathode_plane.h = z_dimension_foils;                     // width and height in cm

  // get scintpoint coords relative to centre of cathode plane
  TVector3 cathodeCentrePoint(x_foils,y_foils,z_foils);
  TVector3 ScintPoint_relative = ScintPoint - cathodeCentrePoint; 

  // calculate solid angle of cathode from the scintillation point
  double solid_angle_cathode = solid(cathode_plane, ScintPoint_relative);

  // calculate distance and angle between ScintPoint and hotspot
  // vast majority of hits in hotspot region directly infront of scintpoint,therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
  double distance_cathode = std::abs(x_foils - ScintPoint[0]); 

  double cosine_cathode = 1;
  double theta_cathode = 0;

  // calculate hits on cathode plane via geometric acceptance
  double cathode_hits_geo = exp(-1.*distance_cathode/L_abs) * (solid_angle_cathode / (4.*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = (theta_cathode/delta_angulo);  
  // correction
  double pars_ini_[4] = {fGHVUVPars[0][j], fGHVUVPars[1][j], fGHVUVPars[2][j], fGHVUVPars[3][j]};
  double GH_correction = GaisserHillas(distance_cathode, pars_ini_);

  double cathode_hits_rec = GH_correction*cathode_hits_geo/cosine_cathode;


  // 2). calculate number of these hits which reach the optical channel from the hotspot via solid angle 
  
  // calculate hotspot location  
  TVector3 v_to_wall(x_foils - ScintPoint[0],0,0);        
  TVector3 hotspot = ScintPoint + v_to_wall;

  // solid angle :
  double solid_angle_detector = 0;
  // rectangular aperture
  if (optical_detector_type == 1) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint; 
    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
    detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window

    
    // get hotspot coordinates relative to detpoint
    TVector3 emission_relative = hotspot - OpDetPoint;

    // calculate solid angle of optical channel
    solid_angle_detector = solid(detPoint, emission_relative);
  }  
  // disk aperture
  else if (optical_detector_type == 0) {
    // offset in z-y plane
    double d = sqrt(pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
    // drift distance (in x)
    double h =  sqrt(pow(hotspot[0] - OpDetPoint[0],2));
    // Solid angle of a disk
    solid_angle_detector = Disk_SolidAngle(d, h, radius);
  }
  else {
    std::cout << "Error: Invalid optical detector type." << endl;
    exit(1);
  }

  // calculate number of hits via geometeric acceptance  
  double hits_geo = (solid_angle_detector / (2*pi)) * cathode_hits_rec;

  /
  // distance to hotspot
  double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
  // distance from hotspot to arapuca
  double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
  // angle between hotspot and arapuca
  double cosine_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
  double theta_vis = acos(cosine_vis)*180./pi;

  // apply correction curves, 5th order polynomial 
  int k = (theta_vis/delta_angle);
  double hits_rec = VIS_pol[k]->Eval(distance_vuv)*hits_geo/cosine_vis;

  // apply border correction
  // calculate radial distance from centre
  double r_distance = sqrt( pow(ScintPoint[1] - y_foils, 2) + pow(ScintPoint[2] - z_foils, 2)); 
  // interpolate in x for each r bin
  std::vector<double> interp_vals = {0,0,0,0};
  for (int i = 0; i < 4; i++){
    interp_vals[i] = interpolate(vDistances_x, VIS_RS60_SP_Borders[k][i], std::abs(ScintPoint[0]), false);
  }
  // interpolate in r
  double border_correction = interpolate(vDistances_r, interp_vals, r_distance, false);
  // apply correction
  double hits_rec_borders = border_correction * hits_rec / cosine_vis;

  // round final result
  int hits_vis = std::round(hits_rec_borders);

  return hits_vis;
}


// gaisser-hillas function definition
Double_t semi_analytic_hits::GaisserHillas(double x,double *par) {
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-x)/par[2]);
  
  return ( Normalization*Term*Exponential);
}


// solid angle of rectanglular aperture calculation functions

double semi_analytic_hits::omega(const double &a, const double &b, const double &d) const{

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double semi_analytic_hits::solid(const acc& out, const TVector3 &v) const{

  //v is the position of the track segment with respect to 
  //the center position of the arapuca window 
 
  // arapuca plane fixed in x direction	

  if( v.Y()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.X());
  }
  
  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  
  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}


// solid angle of circular aperture
double semi_analytic_hits::Disk_SolidAngle(double *x, double *p) {
  const double d = x[0];
  const double h = x[1];
  const double b = p[0];
  if(b <= 0. || d < 0. || h <= 0.) return 0.; 
  const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
  if(d == 0) {
    return 2.*TMath::Pi()*(1.-aa);
  }
  const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
  const double cc = 4*b*d/((b+d)*(b+d));

  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  if(TMath::Abs(ROOT::Math::comp_ellint_1(bb) - bb) < 1e-10 && TMath::Abs(ROOT::Math::comp_ellint_3(cc,bb) - cc) <1e-10) {
    throw(std::runtime_error("please do gSystem->Load(\"libMathMore.so\") before running Disk_SolidAngle for the first time!"));
  }
  if(d < b) {
    return 2.*TMath::Pi() - 2.*aa*(ROOT::Math::comp_ellint_1(bb) + TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb));
  }
  if(d == b) {
    return TMath::Pi() - 2.*aa*ROOT::Math::comp_ellint_1(bb);
  }
  if(d > b) {
    return 2.*aa*(TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb) - ROOT::Math::comp_ellint_1(bb));
  }

  return 0.;
}

double semi_analytic_hits::Disk_SolidAngle(double d, double h, double b) {
  double x[2] = { d, h };
  double p[1] = { b };
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  return Disk_SolidAngle(x,p);
}

double semi_analytic_hits::interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate ) {
  int size = xData.size();
  int i = 0;                                          // find left end of interval for interpolation
  if ( x >= xData[size - 2] )                         // special case: beyond right end
    {
      i = size - 2;
    }
  else
    {
      while ( x > xData[i+1] ) i++;
    }
  double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
  if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
    {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
    }
  double dydx = ( yR - yL ) / ( xR - xL );            // gradient
  return yL + dydx * ( x - xL );                      // linear interpolation
}