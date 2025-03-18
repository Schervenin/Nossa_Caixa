// --------- Spherical Harmonics Sampling

// Author: Mauricio Moralles
// date:   16/11/2020

#include "SphericalHarmonicsSampling.hh"

#include "G4LegendrePolynomial.hh"
#include "Randomize.hh"

#include "globals.hh"

#include <iostream>
#include <iomanip>
//#include <string>
#include <fstream>
//#include <vector>
#include <complex>
#include <cmath>

using namespace std;


// ----- Constructor
SpherHarm::SpherHarm(G4int km)
{
  maxfW = -9.9e15;  // is calculated in the PrepareForSampling method
  minfW =  9.0e15;  // is calculated in the PrepareForSampling method

  kmax = km;
  indmax = kmax+1;  // index from 0 to kmax
  Wkq = new std::complex<G4double>*[indmax];    // k: from 0 to k
  for(G4int i=0; i<indmax; i++)
    Wkq[i]  = new std::complex<G4double>[2*(indmax-1)+1]; // q: -k <= q <= +k
    
  ResetValues();  // initialize with W_00 = (1.0, 0.0) (spherical distrib.)
}

// ----- Destructor
SpherHarm::~SpherHarm() {;}


// ----- Insert Wkq values
void SpherHarm::ReadFromFile(G4String fileName)
{
  std::ifstream WklFile;
  WklFile.open(fileName,std::ios::in);
  if(!WklFile)
  {
      G4cout << "Can not open file " << fileName << endl;
      G4cout << "Only W_00 = (1.0,0.0) will be defined" << endl;
      ResetValues();
      return;
  }
  
  G4int k,l;
  G4double realP, imP;
  G4cout << "\nReading W(k,l) coefficients from file " << fileName << G4endl;
  while( !WklFile.eof() ) 
  {
    WklFile >> k;
    
    if(k<0)
    {
      G4cout << " k < 0. Reading of file " << fileName << " will stop." << endl;
      WklFile.close();
      PrepareForSampling();
      return;
    }
      
    if(k>kmax)
    {
      G4cout << " k > " << kmax << " in file " << fileName << ". File reading will stop." << endl;
      G4cout << " Create a new SpherHarm object with larger kmax if needed. " << endl;
      WklFile.close();
      PrepareForSampling();
      return;      
    }
      
    for(G4int i = 0; i < 2*k+1; i++)
    {        
      if(!(WklFile >> l >> realP >> imP))
      {
        if(WklFile.eof())
        {
          WklFile.close();
          PrepareForSampling();
          return;
        }
        else
        {
          G4cout << " WRONG LINES IN FILE " << fileName << endl;
          G4cout << "   For each k value, 2*k+1 lines must be given" << endl;
          G4cout << "   Only W_00 = (1.0,0.0) will be defined" << endl;
          ResetValues();
          WklFile.close();
          return;
        }
      }
      set(k,l,realP,imP);
//      G4cout << "Read W("<<k<<","<<l<<") = "<<realP<<" , "<<imP<<G4endl;
    }    
  } 
  WklFile.close();
  PrepareForSampling();
  return;
}

// ----- Reset
void SpherHarm::ResetValues()
{
  for(G4int i = 0; i<indmax; i++)
    for(G4int j = 0; j<(2*(indmax-1)+1); j++)
      Wkq[i][j] = std::complex<G4double>(0.0,0.0);
  Wkq[0][0] = std::complex<G4double>(1.0,0.0); // para k=0, m=0
  PrepareForSampling();
  return;
}

// ----- Insere valores
void SpherHarm::set(G4int k, G4int q, G4double re, G4double im)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
  {
    Wkq[k][k+q] = std::complex<G4double>(re,im);
//    G4cout << "k = " << k << " - q = " << q << " - k+q = " << k+q << endl;
    G4cout << " Set Wkq("<<k<<","<<q<<")= " <<Wkq[k][k+q]<<endl;
  }
  else G4cout << "Invalid index: k = " << k << "; q = " << q << endl;
  return;
}

// ----- Lê valores
std::complex<G4double> SpherHarm::get(G4int k, G4int q)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
    return Wkq[k][k+q];
  else
  {
    G4cout << "Invalid index: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

G4double SpherHarm::getReal(G4int k, G4int q)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
    return Wkq[k][k+q].real();
  else
  {
    G4cout << "Invalid index: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

G4double SpherHarm::getIm(G4int k, G4int q)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
    return Wkq[k][k+q].imag();
  else
  {
    G4cout << "Invalid index: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

// Função para gerar os Harmônicos Espéricos
G4double SpherHarm::fSH(G4double th, G4double ph)
{
  G4LegendrePolynomial *LP = new G4LegendrePolynomial();
  G4double PI = 3.14159265359;
  if(th>PI){th=2.0*PI-th; ph=ph+PI;}
  G4double x=cos(th);
  std::complex<G4double> Wc(0.,0.);
//  for(G4int k=0;k<=Kmax;k+=2)
  for(G4int k=0;k<=kmax;k++)    
  {
    for(G4int q=-k;q<=k;q++)
    {
      G4double re,im;
      if(q>=0) {re=cos(q*ph); im=-sin(q*ph);}
      else {re=cos(q*(PI-ph)); im=-sin(q*(PI-ph));}
      std::complex <G4double> cf(re,im);
//      Wc += get(k,q)*gsl_sf_legendre_sphPlm(k,abs(q),x)*cf;
      Wc += get(k,q)*LP->EvalAssocLegendrePoly(k, abs(q), x)*cf;
    }
  }
  return Wc.real();   
}

void SpherHarm::PrepareForSampling()
{
  // Store the maximum and minimum values
  maxfW = -9.9e15;
  minfW =  9.9e15;
  G4double temp  = 0.0;
  G4double PI    = 3.14159265359;
  G4double TwoPI = PI*2.0;
  // the angular step to evaluate maximum and minimum values
  G4double deltaRad = 0.2*PI/180.;  // 0.2 degrees = 0.00349 radians

  for(G4double phi = 0.0+deltaRad/2.0; phi < TwoPI; phi = phi+deltaRad)
  {
    for(G4double theta = 0.0+deltaRad/2.0; theta < PI; theta = theta+deltaRad)
    {
      temp = fSH(theta,phi);
      if(temp > maxfW) maxfW = temp;
      if(temp < minfW) minfW = temp;
    }
  }
//  G4cout << "minfW = " << minfW << "\n";
//  G4cout << "maxfW = " << maxfW << "\n";
  return;  
}

void SpherHarm::SampleThetaPhi(G4double *th, G4double *ph)
{
  G4double PI    = 3.14159265359;
  G4double TwoPI = PI*2.0;
  G4bool success = false;
  while(!success) // sampling using the rejection method
  {
    G4double phiRand   = TwoPI*G4UniformRand();;  // random
    G4double thetaRand = PI*G4UniformRand();;  
    G4double yRand     = maxfW*G4UniformRand();;
    if(yRand < minfW)
    {
      *th = thetaRand;
      *ph = phiRand;
      success = true;
    }
    else if(yRand < fSH(thetaRand,phiRand) )
    {
      *th = thetaRand;
      *ph = phiRand;
      success = true;
    }
  }
  return;
}
// ---------------------------------------------------------

void SpherHarm::printSpherHarmInfo(std::ofstream *output)
{
    *output << G4endl << " --- Spherical Harmonics Coefficients" << G4endl;
    for(G4int i = 0; i<indmax; i++)
      for(G4int j = 0; j<(2*(indmax-1)+1); j++)
        *output << " Wkq("<<i<<","<<j<<")= " <<Wkq[i][j]<<endl;
        
    *output << G4endl;
}

