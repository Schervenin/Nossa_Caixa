// CalorHistoManager.cc

// Description: output data manager for the Calorimeter simulation
// Ntuples are recorded with the ROOT format (.root file)
// The Ntuple contents:  (Event Number, PDG code, Detector Index, Energy, Time)

// Authors: Mauricio Moralles,
// version 4: 02/10/2020 
//            - tree with vectors in one row instead of ntuples with several rows
//            - inclusion of Theta and Phi angles

// version 2: 14/04/2018
  

#include "CalorHistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"

#include "G4AnalysisManager.hh"

using std::vector;

CalorHistoManager::CalorHistoManager()
{
  eventNb = 0;

  // For better perfomance with vectors use "reserve" 
  // Assume that each event will produce less than this number, but it is no problem
  // if it produces more (the vector will increase normally, but at the cost 
  // of longer processing time). 
  G4int nHitsMax = 200;

  PrimPDG.reserve(nHitsMax);
  PDG.reserve(nHitsMax);
  Detector.reserve(nHitsMax);
  Energy.reserve(nHitsMax);
  Time.reserve(nHitsMax);
  Phi.reserve(nHitsMax);
  Theta.reserve(nHitsMax);

  // Clear the vectors
  PrimPDG.clear();
  PDG.clear();
  Detector.clear();
  Energy.clear();
  Time.clear();
  Phi.clear();
  Theta.clear();
}

//______________________________________________________________________________

CalorHistoManager::~CalorHistoManager()
{;}

//______________________________________________________________________________

void CalorHistoManager::Book(G4String name)
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);
  G4cout << "Using analysis type: " << analysisManager->GetType() << G4endl;
//  analysisManager->SetNtupleMerging(true);  // only for root type;
      
  // Open an output file
  //
  fileName = name;
  G4bool fileOpen = analysisManager->OpenFile(fileName);
  if (! fileOpen) {
    G4cerr << "\n---> CalorHistoManager::Book(): cannot open " 
           << analysisManager->GetFileName() << G4endl;
    return;
  }
  
  // Create ntuples.
  // Ntuples ids are generated automatically starting from 0.
  // The start value can be changed by:
  // analysisManager->SetFirstMtupleId(1);  
  
  // Create 1st ntuple (id = 0)
  treeID    = analysisManager->CreateNtuple("tree", "Calorimeter Output");
  eventID   = analysisManager->CreateNtupleIColumn("Event"); // column Id = 0
  NhitsID   = analysisManager->CreateNtupleIColumn("Nhits"); // column Id = 1 
  PrimPDGID = analysisManager->CreateNtupleIColumn("PrimPDG",PrimPDG); // column Id = 2  
  PDGID     = analysisManager->CreateNtupleIColumn("PDG",PDG); // column Id = 3  
  detID     = analysisManager->CreateNtupleIColumn("DetID",Detector); 
  enID      = analysisManager->CreateNtupleDColumn("En",Energy);   
  timeID    = analysisManager->CreateNtupleDColumn("Time",Time);
  phiID     = analysisManager->CreateNtupleDColumn("Phi",Phi); 
  thetaID   = analysisManager->CreateNtupleDColumn("Theta",Theta); 

  analysisManager->FinishNtuple();
  
  G4cout << "\n----> Output file is open in " 
         << analysisManager->GetFileName() << "." 
         << analysisManager->GetFileType() << G4endl;
  
  eventNb = 0;
  
  // Clear the vectors
  PrimPDG.assign(PrimPDG.size(),0);
  PrimPDG.clear();
  PDG.assign(PDG.size(),0);
  PDG.clear();
  Detector.assign(Detector.size(),0);
  Detector.clear();
  Energy.assign(Energy.size(),0);
  Energy.clear();
  Time.assign(Time.size(),0);
  Time.clear();
  Phi.assign(Phi.size(),0);
  Phi.clear();         
  Theta.assign(Theta.size(),0);
  Theta.clear();
}

//______________________________________________________________________________

void CalorHistoManager::Save()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->Write();
  analysisManager->CloseFile(); 
   
  G4cout << "\n----> The tree is saved\n" << G4endl;
      
  return;
}

//______________________________________________________________________________

void CalorHistoManager::FillNtuple()
{
  // FillNtuple is called at the end of each event with Nhits > 0
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  Nhits = PDG.size();
  if(Nhits>0)
  {
    // Not vectorized variables must be filled explicitely
    analysisManager->FillNtupleIColumn(treeID, eventID, eventNb);
    analysisManager->FillNtupleIColumn(treeID, NhitsID, Nhits);
    // The Fill does not apply to vectors. AddNtupleRow fills them automatically
    analysisManager->AddNtupleRow(treeID); // one Row per event, with Nhits PDG, En, Time, Theta, Phi
  }
  // Clear the vectors
  PrimPDG.clear();  
  PDG.clear();
  Detector.clear();
  Energy.clear();
  Time.clear();
  Theta.clear();
  Phi.clear();
  
  Nhits = 0;  // prepare for the enext event
}  

