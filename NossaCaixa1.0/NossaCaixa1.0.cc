
// --------------------------------------------------------------
//      GEANT 4 - Calor.cc
// Program to simulate the Gamma Calirometer of the NUMEN Project
// - Semi-spherical chamber
// - LYSO detectors coupled to SiPMs
// --------------------------------------------------------------
// Author: Mauricio Moralles

// version 5.0 - 25/11/2024
//   - Changes to compile with Geant4 version 11.2

// version 4.0 
//   - inclusion of angular distribution to simulate the demonstrator

// version 3.1 -06/03/2020
//   - Changes in the DetectorElement: inclusion of more details like BaSO4 reflector,
//     Teflon and Aluminum cap with different thicknesses for side and entrance window
//   - DetectorConstruction: now the reference for detector positioning is the 
//     central point of the entrance window. It was the center of the scintillator 
//     in the old versions.

// version 3.0 - 22/06/2019
//   - Inclusion of the option of cylindrical geometry for the detectors
//   - Inclusion of a target holder as described by Calvo in the meeting of April 2019

// version 2.1 - Date 04/05/2018
//   - Inclusion of the option for intrinsic radioactivity in scintillators
//   - Inclusion of a dislocation option to expand the scintillator positions
//     in order to avoid overlaping of volumes

// version 2.0 - Date 16/04/2018
//   Options: Particle, Radioactive Ion, Phase Space file

#include "G4Types.hh"
#include "globals.hh"

#include "G4RunManager.hh"
//#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "CalorDetectorConstruction.hh"
#include "CalorActionInitialization.hh"
#include "CalorPhysicsList.hh"   // from example /advanced/hadrontherapy

// Visual Manager - version 10.06:
// The preprocessor macros G4UI_USE and G4VIS_USE are removed.
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"



int main(int argc,char** argv)
{
  // --- Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc, argv);
  
  // Build the RunManager
  G4RunManager* runManager = new G4RunManager;
//  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // basic classes
  CalorPhysicsList   *calorPhysicsList;  
  CalorDetectorConstruction  *detConstruction;  

  calorPhysicsList           = new CalorPhysicsList();
  runManager->SetUserInitialization(calorPhysicsList);  // before the primary generator
  
  detConstruction          = new CalorDetectorConstruction(); 

  // set initialization classes

  runManager->SetUserInitialization(detConstruction);
  
  auto actionInitialization = new CalorActionInitialization(detConstruction);
  runManager->SetUserInitialization(actionInitialization);

  // ---------------------------------- Initialize visualization
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // Initialize G4 kernel
  // runManager->Initialize(); // this is done in the macro 

  // ------------------- Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // -------------------------------- Process macro or start UI session
  if ( !ui ) 
  { 
    // batch mode; execute an argument macro file if exist
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else 
  { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute initInter.mac");
    ui->SessionStart();
    delete ui;
  }

  // -------------------------------------------- Job termination
  delete visManager;
  delete runManager;

//  return 0;
}

