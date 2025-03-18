
// CalorRunAction.cc
// Author: Mauricio Moralles
// v5: (22/11/2024)
// v1: (28/03/2018)


#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"

#include "CalorRunAction.hh"
#include "CalorRunActionMessenger.hh"
#include "CalorEventAction.hh"
#include "CalorDetectorConstruction.hh"
#include "CalorHistoManager.hh"
#include "CalorPrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"

#include "globals.hh"

#include <fstream>  // for file management

//_________________________________________________________________________________
CalorRunAction::CalorRunAction(CalorDetectorConstruction* DC, CalorHistoManager* HM, CalorEventAction * EA)
:G4UserRunAction(),
DetCon(DC),HistMan(HM),EventAct(EA)
{
  G4cout << "Building the CalorRunAction..." << G4endl;
  outputFileName = "Calor_01.root";
  
  RunActMessenger = new CalorRunActionMessenger(this);
}

//_________________________________________________________________________________
CalorRunAction::~CalorRunAction(){
  delete RunActMessenger;
}


//_________________________________________________________________________________
void CalorRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "# Run " <<  aRun->GetRunID() << " start." << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // An instance of the AnalysisManager has been already initialized at the HistMan
  
  HistMan->Book(outputFileName);
  
  G4int eventsToBeProcessed = aRun->GetNumberOfEventToBeProcessed();
  EventAct->SetTotalRunEvtsPtr(eventsToBeProcessed);
  G4cout << " Number of Events for this Run: " << eventsToBeProcessed << G4endl; 
  G4RunManager *RM = G4RunManager::GetRunManager();
  CalorPrimaryGeneratorAction *PG = (CalorPrimaryGeneratorAction*)RM->GetUserPrimaryGeneratorAction();
  PG->GetNumberOfEvents(eventsToBeProcessed);
}

//_________________________________________________________________________________
void CalorRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

  // Get the beam particle and energy
  G4RunManager *RM = G4RunManager::GetRunManager();
  CalorPrimaryGeneratorAction *PG = (CalorPrimaryGeneratorAction*)RM->GetUserPrimaryGeneratorAction();
    
  G4cout << " End of Run                  " << aRun->GetRunID() << G4endl;
  G4cout << " Number of Events:           " << NbOfEvents << G4endl;
  if(PG->GetPrimaryTypeCode() == 3)
  {
    G4cout << " Primary particles from Phase Space File." << G4endl;
    G4cout << " Phase Space Directory: " << PG->GetPhaseSpaceFileDirectory() << G4endl;
    G4cout << " Phase Space File Name: " << PG->GetPhaseSpaceFileName() << G4endl;
    if(PG->GetEndOfPhaseSpaceFile())
    {
      G4cout<< " End of Phase Space File reached." << G4endl;
      G4cout << " All events of the Phase Space File were used." << G4endl;
    }    
  }
    
  G4cout << " Saving results to root file " << analysisManager->GetFileName() << G4endl;
  HistMan->Save();   

  // Description of this run saved in text file
  G4cout << "--> Saving Run info to file CalorRunInfo.txt" << G4endl << G4endl;

  std::ofstream outFile;
  outFile.open("CalorRunInfo.txt",std::ios::out);
  outFile << " End of Run           " << aRun->GetRunID() << G4endl;
  outFile << " Number of Events:    " << NbOfEvents << G4endl;
  outFile << " Data file:           " << outputFileName << G4endl << G4endl;
    
  DetCon->printDetectorInfo(&outFile);
  PG->printPrimaryGeneratorInfo(&outFile);
  // Extra info in the case of intrinsic activity simulation
  G4int       nScintillators;
  G4double    detVolume, calorScintiVolume, calorActivity;
  G4Material* matScint = DetCon->GetDetectorElement()->getMaterialScintillator();
  G4double    LYSO_activity  = 307.0;  // Bq/cm3
  G4double    LaBr3_activity = 1.5;    // Bq/cm3
  G4double    timeEquivalence;
  if(PG->GetPrimaryTypeCode() == 4)  // Intrinsic Radiation
  {
    nScintillators    = DetCon->GetNumberOfDetectors();
    G4String shape    = DetCon->GetDetectorShape();
    G4double height   = DetCon->GetScintillatorHeight();
    G4double section  = DetCon->GetScintillatorSection();
    
    // Volume of one scintillator
    if(shape == "Box") detVolume = section*section*height;
    else detVolume = M_PI*pow(section/2.0,2)*height;
    // Total volume of all scintillators
    calorScintiVolume = (G4double)(nScintillators)*detVolume;
       
    // Total activity of the calorimeter
    if(matScint->GetName() == "LYSO")
      calorActivity     = (calorScintiVolume/cm3)*LYSO_activity;
    else if(matScint->GetName() == "LaBr3")
      calorActivity     = (calorScintiVolume/cm3)*LaBr3_activity;
    else
      calorActivity     = 0.0;
       
    // Time corresponding to the chosen number of primaries 
    // i.e. number of decays of Lu-176 (LYSO) or La-138 (LaBr3)
    timeEquivalence = (G4double)(NbOfEvents)/calorActivity;
       
    // ----- Print on the screen
    G4cout << "The calorimeter has     " << nScintillators << " of " <<  matScint->GetName() << G4endl;
    G4cout << "Each scintillator has   " << detVolume/cm3 << " cm3" << G4endl; 
    G4cout << "The calorimeter has     " << calorScintiVolume/cm3 << " cm3 of scintillator" << G4endl;
    if(matScint->GetName() == "LYSO")
      G4cout << "LYSO has activity of    " << LYSO_activity << " Bq/cm3" << G4endl;
    else if(matScint->GetName() == "LaBr3")
      G4cout << "LaBr3 has activity of   " << LaBr3_activity << " Bq/cm3" << G4endl;
    G4cout << "Calorimeter activity:   " << calorActivity << " Bq" << G4endl;
    G4cout << "Number of primaries:    " << NbOfEvents << " corresponds to " << timeEquivalence << " seconds of measurement" << G4endl;
       
    // ----- Print in the report file
    outFile << "The calorimeter has   " << nScintillators << " of " <<  matScint->GetName() << G4endl;
    outFile << "Each scintillator has " << detVolume/cm3 << " cm3" << G4endl; 
    outFile << "The calorimeter has   " << calorScintiVolume/cm3 << " cm3 of scintillator" << G4endl;
    if(matScint->GetName() == "LYSO")
      outFile << "LYSO has activity of  " << LYSO_activity << " Bq/cm3" << G4endl;
    else if(matScint->GetName() == "LaBr3")
      outFile << "LaBr3 has activity of  " << LaBr3_activity << " Bq/cm3" << G4endl;
    outFile << "Calorimeter activity: " << calorActivity << " Bq" << G4endl;
    outFile << "Number of primaries:  " << NbOfEvents << " corresponds to " << timeEquivalence << " seconds of measurement" << G4endl;
  }
  outFile.close();
  
  HistMan->ResetEventNumber();
  
  
  return;
}

//_________________________________________________________________________________
void  CalorRunAction::SetOutputFileName(G4String name) 
{
    outputFileName = name;
}

//_________________________________________________________________________________

