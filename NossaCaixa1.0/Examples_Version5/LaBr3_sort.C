// --------------- Function to sort te energy values according
// to the LaBr3 resolution
// The FWHM function is FWHM = a1 + a2*sqrt(Energy)
// parameters a1 and a2 depends on the detector

// use: .L LaBr3_sort.C


// 511 keV:  FWHM = 4 channels 
// 14 channels from 480 keV to 550 keV 
// ~ 70 keV/14 ch = 5 keV/ch -> FWHM = 4*5 = 20 keV
// -> resolution ~ 20/511 = 3.9%

// 1275 keV: FWHM = 6.5 channels
// 18 channels from 1220 keV to 1310 keV
// ~ 90 keV/18 channels ~ 5 keV/ch -> FWHM = 6.5*5 = 32.5 keV
// -> resolution = 32.5/1275 = 2.5%

// Since we have only two points and one is for the 511 keV,
// wich is always wider than standard lines, we used the
// values of a Saint-Gobain detector corrected for the difference
// presented at 1275 keV

//a1 = 0.992351
//a2 = 0.962132



TRandom ran;

Double_t LaBr3_sort(Double_t Energy)
{
  Double_t a1 = 1.84051;
  Double_t a2 = 1.01567;
  Double_t sigma = (a1+a2*sqrt(Energy))/2.35;  // FWHM = 2.35*sigma
    
  return ran.Gaus(Energy,sigma);
};
 
