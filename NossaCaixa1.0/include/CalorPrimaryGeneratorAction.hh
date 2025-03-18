// CalorPrimaryGeneratorAction.hh

// Author: Mauricio Moralles
// Date:   16/11/2020

// 16/11/2020: Inclusion of angular distribution for the Particle type

#ifndef CalorPrimaryGeneratorAction_h
#define CalorPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"

#include <vector>

class CalorPrimaryGeneratorMessenger;
class G4ParticleGun;
class G4Event;
class G4ParticleDefinition;
class G4RootAnalysisReader;
class G4ParticleTable;
class CalorDetectorConstruction;
class SpherHarm;


class CalorPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
      CalorPrimaryGeneratorAction(CalorDetectorConstruction* DC);
      ~CalorPrimaryGeneratorAction() override;

  public:
      void GeneratePrimaries(G4Event* anEvent);

  private:
      G4int                   typeOfPrimary;    // 1:Gamma; 2: RadioIsotope; 3: Phase Space File
      G4String                StrTypeOfPrimary; // 1:Gamma; 2: RadioIsotope; 3: Phase Space File
      
      // General
      G4ParticleGun*          particleGun = nullptr;
      G4ParticleTable*        particleTable = nullptr;
      G4ParticleDefinition*   primaryParticleDef = nullptr;
      G4String                particleName;
      G4ThreeVector           particlePosition;
      G4ThreeVector           particleDirection;
      G4double                particleEnergy;
      G4double                particleTime;
      G4double                phiMin, phiMax, thetaMin, thetaMax;
      G4bool                  angDistFlag;  // to turn on/off the angular distribution sampling
      G4String                WkqDirectory, WkqName, WkqFile;
      G4int                   kMax;
      SpherHarm*              AngDist  = nullptr;
      
      // For phase space file
      G4RootAnalysisReader*   analysisReader = nullptr;
      G4bool                  analysisReaderInstanceFlag;
      G4String                phaseSpaceDirectory, phaseSpaceName, phaseSpaceFile;
      G4bool                  endOfPhaseSpaceFileFlag;
      G4int                   ntupleID;
      G4double                dirX, dirY, dirZ, time;
      G4double                radiusOfThePhaseSpaceSphere;
      G4int                   numberOfEvents;  // for PhaseSpace events
      G4int                   reactionNumber, lastReactionNumber, Nhits;
      G4int                   PDGencoding, volumeNumber;
      std::vector<G4int>      PDGcodeVector, volumeNbVector;      
      std::vector<G4double>   enVector,dXVector,dYVector,dZVector,timeVector;
      
      // For Radioisotope
      G4int                   Z_number, A_number;
      G4double                excitationEnergy;
      
      // For detector intrinsic radioactivity
      G4int                   detShape;  // 1 = Box; 2 = Cylinder
      G4ThreeVector           scintPosition;
      G4double                scintDimX, scintDimY, scintDimZ, scintRadius;
      

  public:
      G4ParticleGun* GetParticleGun()   {return particleGun;}
      
      void           SetPhiMin(G4double pMin)       {phiMin = pMin;}
      void           SetPhiMax(G4double pMax)       {phiMax = pMax;}
      void           SetThetaMin(G4double tMin)     {thetaMin = tMin;}
      void           SetThetaMax(G4double tMax)     {thetaMax = tMax;}
      
      void           SetPrimaryType(G4int type);
      G4int          GetPrimaryTypeCode()           {return typeOfPrimary;}
      G4String       GetPrimaryType()               {return StrTypeOfPrimary;}
      void           GetNumberOfEvents(G4int ne)    {numberOfEvents = ne;}
      G4int          GetReactionNumber()            {return reactionNumber-1;}
      G4int          GetVolumeNumber()              {return volumeNumber;}
      G4int          GetPrimaryPDGencoding()        {return PDGencoding;}
      G4double       GetPrimaryTime()               {return particleTime;}
      
      void      SetAngDistFlag(G4bool flag);
      void      SetKmax(G4int max)                {kMax = max;}      
      G4bool    GetAngDistFlag()                  {return angDistFlag;}
      void      SetAngDistDirectory(G4String dir);
      G4String  GetAngDistDirectory()             {return WkqDirectory;}
      void      SetAngDistName(G4String name);
      G4String  GetAngDistName()                  {return WkqName;}
      
      
      G4String       GetPhaseSpaceFileDirectory()   {return phaseSpaceDirectory;}
      G4String       GetPhaseSpaceFileName()        {return phaseSpaceName;}
      G4bool         GetEndOfPhaseSpaceFile()       {return endOfPhaseSpaceFileFlag;}
      void           SetRadiusOfSphere(G4double rd) {radiusOfThePhaseSpaceSphere = rd;}
      void           SetPhaseSpaceFileDirectory(G4String dir);
      void           SetPhaseSpaceName(G4String name);
      
      void           SetIon(G4int Z,G4int A, G4double excitationEnergy);
      // to print info at the end of the Run
      void   printPrimaryGeneratorInfo(std::ofstream *);
     
  private:
      CalorPrimaryGeneratorMessenger* primaryMessenger = nullptr;
      CalorDetectorConstruction*      detCon = nullptr;

};

#endif
