

// --------- Spherical Harmonics Sampling

// Author: Mauricio Moralles
// date:   16/11/2020
#include "globals.hh"

#include <complex>

class ofstream;

class SpherHarm
{
 public:
   SpherHarm(G4int km);
   ~SpherHarm();

   void                   set(G4int, G4int, G4double, G4double);
   G4int                  getKmax() {return kmax;}
   std::complex<G4double> get(G4int,G4int);
   G4double               getReal(G4int,G4int);
   G4double               getIm(G4int,G4int);
   
   G4double  fSH(G4double th, G4double ph);
   void      ReadFromFile(G4String fileName);
   void      ResetValues();
   void      PrepareForSampling();
   void      SampleThetaPhi(G4double *th, G4double *ph);
   void      printSpherHarmInfo(std::ofstream *);
  
 private:
   G4int                    kmax, indmax;
   std::complex<G4double> **Wkq;
   G4double                 maxfW, minfW;        // for sampling
};
 
