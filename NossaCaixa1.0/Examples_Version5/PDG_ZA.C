// Utility for PDG output

// Use: .L PDG_ZA.C

// PDG for ions
// 1000ZZAAAI ; I: excited state

//  Int_t ZZZ     = Int_t(ZZZAAAI/10000);
//  Int_t AAAI    = Int_t(ZZZAAAI-ZZZ*10000);
//  Int_t AAA     = Int_t(ZZZAAAI-ZZZ*10000)/10;
//  Int_t AAA     = Int_t(ZZZAAAI-ZZZ*10000)/10;  

#include "TString.h"

Int_t PDG_ZA(Int_t PDG)  // return ZZAAA
{
  if(PDG<1000000000)
    return PDG;
 // nuclei PDG: 1000ZZAAAE
  else
    return PDG-1000000000;
}

Int_t PDG_Z(Int_t PDG)   // return Z
{
  if(PDG<1000) return (-1000-abs(PDG)); // not hadron
  else if(PDG == 2112) return 0;        // neutron
  else if(PDG == 2212) return 1;        // proton, positron
  else
    return Int_t((PDG-1000000000)/10000);
}

Int_t PDG_A(Int_t PDG)   // return A
{
  if(PDG<1000) return (-1000-abs(PDG)); // not hadron
  else if(PDG == 2112) return 1;        // neutron
  else if(PDG == 2212) return 1;        // proton
  else 
  {
    Int_t ZZZAAAI = (PDG-1000000000);
    return  Int_t(ZZZAAAI-(Int_t(ZZZAAAI/10000))*10000)/10;
  }
}
