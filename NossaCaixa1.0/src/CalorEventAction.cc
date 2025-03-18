
// CalorEventAction.cc


#include "CalorEventAction.hh"
#include "CalorHistoManager.hh"
#include "CalorPrimaryGeneratorAction.hh"
#include "CalorSteppingAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"   // to take a pointer to RunManager
#include "G4SystemOfUnits.hh"

#include <vector>

//______________________________________________________________________

CalorEventAction::CalorEventAction(CalorHistoManager* HM, CalorPrimaryGeneratorAction* PG)
:G4UserEventAction(),
 HistMan(HM),PrimGen(PG)
{
  eventNumber = HM->GetEventNb();
  fPrintModulo = 10;
}

//______________________________________________________________________

//CalorEventAction::~CalorEventAction() { }

//______________________________________________________________________

void CalorEventAction::BeginOfEventAction(const G4Event*)
{
  // Begin of Event. The SteppingAction reference variables must be reset
  // Get a pointer to the SteppingAction
  G4RunManager *RM = G4RunManager::GetRunManager();
  CalorSteppingAction *SA = (CalorSteppingAction*)RM->GetUserSteppingAction();
  SA->ResetForNewEvent();

  // Print events progress in 10% units
  if(*eventNumber%fPrintModulo == 0) // fPrintModulo = totalRunEvents/10
    G4cout << "---> Beginning event: " << *eventNumber+1 << " of " << totalRunEvents <<  G4endl;
  return;
}

//______________________________________________________________________

void CalorEventAction::EndOfEventAction(const G4Event* evt)
{
//  if(PrimGen->GetPrimaryTypeCode() == 3) // if primaries from phase space 
//    *eventNumber = PrimGen->GetReactionNumber();
//  else *eventNumber = evt->GetEventID(); 
  
/*  
// TEST: print some information
G4cout << "\n Energy of this event: " << G4endl;
std::vector<G4double>* Ener=HistMan->GetEnergy();
std::vector<G4int>*    PDG=HistMan->GetPDG();
std::vector<G4int>*    Prim=HistMan->GetPrimPDG();
G4int                  Nh=Ener->size();

G4cout << " Number of hits: " << Nh << G4endl;
for(G4int i=0; i< Nh;i++)
  G4cout << " hit: " << i 
  << "  Energy: " << Ener->at(i)/keV << " keV" 
  << " Primary: " << Prim->at(i) 
  << " PDG: " << PDG->at(i) << G4endl;
// end TEST ---------
*/

  *eventNumber = evt->GetEventID();      
  // Fill the ntuples with all entries of this event
  HistMan->FillNtuple();
  return;
}  

//______________________________________________________________________
