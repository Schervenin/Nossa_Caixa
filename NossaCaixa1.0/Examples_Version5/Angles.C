// Convert versor directions to angles

// Use: .L Angles.C

#include "TMath.h"

// TMath::Pi()      = 3.1415926536
// 2*TMath::Pi()    = 6.2831853072
// 0.5*TMath::Pi()  = 1.5707963268
// 1.5*TMath::Pi()  = 4.7123889804
// 180./TMath::Pi() = 57.295779513

Double_t ThetaRad(Double_t dZ)  // radians
{
  return acos(dZ);
}

Double_t ThetaDeg(Double_t dZ)  // degrees
{
  return acos(dZ)*57.295779513;
}

Double_t PhiRad(Double_t dX,Double_t dY)  // radians
{
  if(dX<0)      return atan(dY/dX) + 3.1415926536;  // 2nd and 3rd quadrants
  else if(dY>0) return atan(dY/dX);                 // 1st quadrant
  else if(dY<0) return atan(dY/dX) + 6.2831853072;  // 4th quadrant
  else if(dX==0.0 && dY==1.0)  return 1.5707963268; // pi/2
  else if(dX==0.0 && dY==-1.0) return 4.7123889804; // 3*pi/2
  else return -9999.9;  // problem with input values
}

Double_t PhiDeg(Double_t dX,Double_t dY)  // degrees
{
  if(dX<0)      return atan(dY/dX)*57.295779513 + 180.; // 2nd and 3rd quadrants
  else if(dY>0) return atan(dY/dX)*57.295779513;        // 1st quadrant
  else if(dY<0) return atan(dY/dX)*57.295779513 + 360.; // 4th quadrant
  else if(dX==0.0 && dY==1.0)  return 90.;              // pi/2
  else if(dX==0.0 && dY==-1.0) return 270.;             // 3*pi/2
  else return -9999.9;  // problem with input values
}





