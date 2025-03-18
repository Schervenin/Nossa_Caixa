// CalorHistoManager.hh

// Description: output data manager for the Calorimeter simulation
// Ntuples are recorded with the ROOT format (.root file)
// The Ntuple contents:  (Event Number, Detector Index, Energy, Time)

// Authors: Mauricio Moralles,
// v4.0: (02/10/2020): now the Ntuple are formed by rows with vectors (tree)  

// v1: (28/03/2018)

#ifndef CalorHistoManager_h
#define CalorHistoManager_h 1

#include "globals.hh"

#include <vector>
//#include "g4root.hh"   // deprecated after Version 11
//#include "g4csv.hh"
//#include "g4xml.hh"

//______________________________________________________________________________

class CalorHistoManager
{
  public:
    CalorHistoManager();
   ~CalorHistoManager();

    void Book(G4String name);
    void Save();
    

    void FillNtuple();
    
    void SetFileName(G4String name)   {fileName = name;}
    
  private:
    // Name of output file
    G4String fileName;
    // Ntuple ID
    G4int treeID;
    
    G4int eventNb, Nhits;   // number of event, hits of each event
    // Columns IDs
    G4int eventID, NhitsID, PrimPDGID, PDGID, detID, enID, timeID;
    G4int thetaID, phiID;
    
    std::vector<G4int>    PrimPDG;
    std::vector<G4int>    PDG;
    std::vector<G4int>    Detector;
    std::vector<G4double> Energy;
    std::vector<G4double> Time;
    std::vector<G4double> Theta;    
    std::vector<G4double> Phi;
    
public:
    G4int*                 GetEventNb()  { return &eventNb;  } 
    std::vector<G4int>*    GetPrimPDG()  { return &PrimPDG;  }
    std::vector<G4int>*    GetPDG()      { return &PDG;      }
    std::vector<G4int>*    GetDetector() { return &Detector; }
    std::vector<G4double>* GetEnergy()   { return &Energy;   }
    std::vector<G4double>* GetTime()     { return &Time;     }
    std::vector<G4double>* GetTheta()    { return &Theta;    }    
    std::vector<G4double>* GetPhi()      { return &Phi;      }
    
    void  ResetEventNumber() { eventNb = 0; }
};

//______________________________________________________________________________

#endif

