#ifndef TIMEPARAMSERVICES_H
#define TIMEPARAMSERVICES_H

// ======================================================================
// Contains various functions required in generaing the timing parameterisations in the timeparam class

//   Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)

#include "TF1.h"
#include <vector>



double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate )
{
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

double* interpolate( const std::vector<double> &xData, const std::vector<double> &yData1, const std::vector<double> &yData2,
		  const std::vector<double> &yData3, double x, bool extrapolate)
{
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
  double xL = xData[i], xR = xData[i+1];// points on either side (unless beyond ends)
  double yL1 = yData1[i], yR1 = yData1[i+1], yL2 = yData2[i], yR2 = yData2[i+1], yL3 = yData3[i], yR3 = yData3[i+1]; 
  
  if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
    {
      if ( x < xL ) {yR1 = yL1; yR2 = yL2; yR3 = yL3;}
      if ( x > xR ) {yL1 = yR1; yL2 = yR2; yL3 = yR3;}
    }
  double dydx1 = ( yR1 - yL1 ) / ( xR - xL );            // gradient
  double dydx2 = ( yR2 - yL2 ) / ( xR - xL );
  double dydx3 = ( yR3 - yL3 ) / ( xR - xL );

  double *yy = new double[3]; 
  yy[0] = yL1 + dydx1 * ( x - xL );// linear interpolations
  yy[1] = yL2 + dydx2 * ( x - xL );
  yy[2] = yL3 + dydx3 * ( x - xL );
  
  return yy;                      
}

/*                       ======TIMING PARAMETRIZATION=====               */
double finter_d(double *x, double *par) {
  
  double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
  double y2 = TMath::Exp(par[3]+x[0]*par[4]);
  
  return TMath::Abs(y1 - y2);
}

Double_t model_close(Double_t *x, Double_t *par)
{
  // par0 = joining point
  // par1 = Landau MPV
  // par2 = Landau width
  // par3 = normalization
  // par4 = Expo cte
  // par5 = Expo tau
  // par6 = t_min
  
  Double_t y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
  Double_t y2 = TMath::Exp(par[4]+x[0]*par[5]);
  if(x[0] <= par[6] || x[0] > par[0]) y1 = 0.;
  if(x[0] < (par[0])) y2 = 0.;
  
  return (y1 + y2);
}

Double_t model_far(Double_t *x, Double_t *par)
{
  // par1 = Landau MPV
  // par2 = Landau width
  // par3 = normalization
  // par0 = t_min

  Double_t y = par[3]*TMath::Landau(x[0],par[1],par[2]);
  if(x[0] <= par[0]) y = 0.;
  
  return y;
}

#endif
