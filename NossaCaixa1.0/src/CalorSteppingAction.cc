// CalorSteppingAction.cc


#include "CalorSteppingAction.hh"

#include "CalorDetectorConstruction.hh"
#include "CalorEventAction.hh"
#include "CalorPrimaryGeneratorAction.hh"
#include "CalorHistoManager.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

#include <vector>

using namespace std;

//_________________________________________________________________________

CalorSteppingAction::CalorSteppingAction(CalorDetectorConstruction* DC, CalorPrimaryGeneratorAction* PG, CalorEventAction* EA, CalorHistoManager* HM)
: G4UserSteppingAction(), 
  DetCon(DC), PrimGen(PG), EventAc(EA), HistMan(HM)
{ 
  scintillatorName = DC->GetScintillatorName();
    
  // set the pointers to the HistoManager vectors
  PrimPDG   = HistMan->GetPrimPDG();
  PDG       = HistMan->GetPDG();
  Detector  = HistMan->GetDetector();
  Energy    = HistMan->GetEnergy();
  Time      = HistMan->GetTime();
  Theta     = HistMan->GetTheta();
  Phi       = HistMan->GetPhi();

  G4cout << "  Got vectors to store the data" << G4endl;    
}

//_________________________________________________________________________

//CalorSteppingAction::~CalorSteppingAction(){ }

//_________________________________________________________________________

void CalorSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  
  G4double depEnergy = aStep->GetTotalEnergyDeposit();
  G4Track* theTrack    = aStep->GetTrack();
  // The primary of all subsequent secondary particles is the one
  // with ParentID = 0
  if(theTrack->GetParentID() == 0)
    primaryPDG = theTrack->GetParticleDefinition()->GetPDGEncoding();
    
// TEST:
//G4cout << "TrackID:  " << theTrack->GetTrackID() << G4endl;
//G4cout << "ParentID: " << theTrack->GetParentID() << G4endl;
  
  if(depEnergy==0.0) return;
  
  // get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
//  G4String volumeName = volume->GetName();
  G4String LogicalVolumeName = volume->GetLogicalVolume()->GetName();
  
  
  // The names of all detectors begin with "D_", defined in the DetectorConstruction
  // The name of the scintillators begin with "D_Scinti". In the next
  // line we compare the first 7 characters of the volume name for this Step
  // and the name of the LogicalVolume of the detectors. If it matches with
  // the scintillator's name, the energy is registered
  // No other volume will have a name beginning with "D_Scint"
  // OBS: in the old version, the detector was composed by only a scintillator, so we coud take the PhysicalVolume name. In this version, 
  // the detector is composed of several LogicalVolumes and the scintillator is only one of them.
    
  // When the Step is inside any of the calorimeter scintillators:
  if (LogicalVolumeName.compare(0,7,scintillatorName,0,7) == 0)
  {
    // get the index of the detector where the energy was deposited
    G4int detectorNumber = aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1);
        
    //G4Track* theTrack    = aStep->GetTrack();
    G4int    PDGencoding = theTrack->GetParticleDefinition()->GetPDGEncoding();
    
  // TEST: print some information

//G4cout << "LogicalVolumeName: " << LogicalVolumeName << G4endl;
//G4cout << "Name: " << LogicalVolumeName << " - Number: " << detectorNumber << G4endl;
//G4cout << "Energy this hit:   " << depEnergy << G4endl;
    
//G4int PrimaryPDG = PrimGen->GetPrimaryPDGencoding();
    // if the type is PhaseSpace, we need to get the particle PDG from the ROOT PDG

//G4cout << "PrimGen PDG:       " << PrimaryPDG << G4endl;  // TEST
//G4cout << "primary PDG:       " << primaryPDG << G4endl;  // TEST
//G4cout << "PDG of this hit:   " << PDGencoding << G4endl; // TEST   
    
  // If electron (PDG = 11) or positron (PDG = -11) as secondary of gamma, attributes energy to PDG of gamma (22)
    if((PDGencoding == 11 || PDGencoding == -11) && theTrack->GetParentID() != 0) 
        PDGencoding = 22;  // gamma
        
//G4cout << "PDG after filter for gamma:  " << PDGencoding << G4endl;     // TEST       
        
    // if this detector was used in the last step (most frequent)
    // the only action to do is to add the energy of this step
    if( (detectorNumber == lastDetectorHit) &&  (PDGencoding == lastPDGHit) )
    {
      Energy->at(indexOfLastHit) += depEnergy;  // adds the energy of this step
//G4cout << "Stepping: add energy to detector " << lastDetectorHit  << " indexed as " << indexOfLastHit  << G4endl;  // for test
      return;
    }
    else  // if not the last detector or not the last PDG, look for other already fired detector
    {
      for(unsigned int index = 0; index < Detector->size(); index++)
      {
        if( (detectorNumber == Detector->at(index)) && (PDGencoding == PDG->at(index)) )
        {
          Energy->at(index) += depEnergy;     // adds the energy of this step
          lastDetectorHit  = detectorNumber; // updates last value
          lastPDGHit       = PDGencoding;    // updates last value
          indexOfLastHit   = index;          // updates last value
//G4cout << "Stepping: add energy to detector " << Detector->at(index)  << " indexed as " << index  << G4endl;  // for test
          return;
        }
      } // if not satisfied above, it is a new pair (detectorNumber,PDGencoding)
      // The registered time, tehta and phi for each pair (detector,PDG) correspond to that of the first occurrence. They are not updated in the next occurrences of the pair.
      // NOTE: Theta and Phi need some algorithm of energy weighted calculation in the case of a position sensitive detector
      
      // TEST
//      G4cout << "PDGs that go to the ROOT file: " << G4endl;
//      G4cout << " PrimaryPDG: " << PrimaryPDG << "  detector PDG: " << PDGencoding << G4endl;
//      G4cout << "PDG that do not go to the ROOT file: " << G4endl;
//      G4cout << " primaryPDG: " << primaryPDG << G4endl;
      
      G4double detectorTime = theTrack->GetGlobalTime();            
      Detector->push_back(detectorNumber);
      PrimPDG->push_back(primaryPDG);
      PDG->push_back(PDGencoding);
      Energy->push_back(depEnergy);
      Time->push_back(detectorTime);
      // Phi and Theta
      G4ThreeVector position = aStep->GetPreStepPoint()->GetPosition();
      G4double x = position.x();
      G4double y = position.y();
      G4double z = position.z();
      if(x==0.0)
      {
        if(y>0.0) Phi->push_back( 0.5*M_PI/degree);
        else      Phi->push_back(1.5*M_PI/degree);
      }
      else 
      {
        G4double PhiTemp = atan(y/x)/degree;
        if(PhiTemp<0.0) Phi->push_back(180.0+PhiTemp);
        else            Phi->push_back(PhiTemp);
      }
      Theta->push_back(acos(z/sqrt(x*x+y*y+z*z))/degree);

      lastDetectorHit = detectorNumber;  // updates last value
      lastPDGHit      = PDGencoding;     // updates last value
      indexOfLastHit  = Detector->size()-1;
      return;
    }
  }
  return;
}

//_________________________________________________________________________

void CalorSteppingAction::ResetForNewEvent()
{
  // values that these variables will never receive during StepAction 
  lastDetectorHit = -9999999;
  lastPDGHit      = -9999999;
  indexOfLastHit  = -9999999;
}
  



