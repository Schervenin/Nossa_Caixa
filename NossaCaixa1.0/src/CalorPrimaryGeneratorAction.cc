
// CalorPrimaryGeneratorAction.cc


#include "CalorPrimaryGeneratorAction.hh"
#include "CalorPrimaryGeneratorMessenger.hh"
#include "CalorDetectorConstruction.hh"
#include "DetectorElement.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4RootAnalysisReader.hh"

//#include "g4root.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericIon.hh"

#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include "SphericalHarmonicsSampling.hh"

#include "globals.hh"

#include <cmath> // to use math functions

// Geant4 version 11 way to use analysis reader
#include "G4RootAnalysisReader.hh"
using G4AnalysisReader = G4RootAnalysisReader;

CalorPrimaryGeneratorAction::CalorPrimaryGeneratorAction(CalorDetectorConstruction* DC)
:detCon(DC)
//particleGun(0),particleTable(0),primaryParticleDef(0),
//AngDist(0),
//analysisReader(0),primaryMessenger(0)
{
  G4cout << G4endl << "Building the PrimaryGeneratorAction... " << G4endl;
    
  particleGun = new G4ParticleGun();
  phiMin   = 0.;
  phiMax   = 2.*M_PI;
  thetaMin = 0.;
  thetaMax = M_PI;
    
  // Initial parameters. Changes can be done via macro file
  particleGun->SetNumberOfParticles(1);
  particleName       = "gamma";
  particleTable      = G4ParticleTable::GetParticleTable();
  primaryParticleDef = particleTable->FindParticle(particleName);
  particlePosition   = G4ThreeVector(0.,0.,0.);
  particleDirection  = G4ThreeVector(0.,0.,1.);
  particleEnergy     = 1.0*MeV;
  particleTime       = 0.;
  
  // Angular distribution
  angDistFlag        = false;  // isotropic by default 
  kMax               = 4;      // maximum k 
  AngDist            = new SpherHarm(kMax);
  WkqDirectory = "./";
  WkqName      = "Wkq_example.dat"; 
  WkqFile      = WkqDirectory+WkqName;
    
  particleGun->SetParticleDefinition(primaryParticleDef);
  particleGun->SetParticlePosition(particlePosition);
  particleGun->SetParticleMomentumDirection(particleDirection);
  particleGun->SetParticleEnergy(particleEnergy);

  // Initial values that can be used in case of changing the typeOfPrimary
  phaseSpaceDirectory = "./";
  phaseSpaceName      = "PhSpace_01.root"; 
  phaseSpaceFile      = phaseSpaceDirectory+phaseSpaceName;
    
  // The next value must be taken from the PhaseSpace simulation: the inner radius 
  // defined in the PhSpDetectorConstruction
  radiusOfThePhaseSpaceSphere = 10.0*cm;
    
  // For radioisotope, typeOfPrimary = 2
  // If you put some kinetic energy to Ba-137 at an excited state
  // it is possible to see the Doppler shift of the 662 keV gamma
  Z_number = 55; A_number = 137;   // Cs-137
  excitationEnergy = 0.;       // ground state

//    Z_number = 56; A_number = 137;   // Ba-137, daughter of Cs-137
//    excitationEnergy = 663.*keV;     // excited state

//    Z_number = 27; A_number = 60;    // Co-60
//    excitationEnergy = 0;     // excited state


  // choice of the default primary type
  typeOfPrimary = 1; // 1:Particle; 2: Isotope; 3: Phase Space File; 4: Intrinsic Radioac.

  analysisReaderInstanceFlag = 0;
  endOfPhaseSpaceFileFlag = 0;
  detShape = 1;  // default value for shape: Box = 1; Cylinder = 2
    
  SetPrimaryType(typeOfPrimary);

  // For better perfomance with vectors use "reserve" 
  // Assume that each event will produce less than this number, but it is no problem
  // if it produces more (the vector will increase normally, but at the cost 
  // of longer processing time). 
  G4int nHitsMax = 200;

  PDGcodeVector.reserve(nHitsMax);
  volumeNbVector.reserve(nHitsMax);
  enVector.reserve(nHitsMax);
  dXVector.reserve(nHitsMax);
  dYVector.reserve(nHitsMax);
  dZVector.reserve(nHitsMax);
  timeVector.reserve(nHitsMax);

  // Clear the vectors
  PDGcodeVector.clear();
  volumeNbVector.clear();
  enVector.clear();
  dXVector.clear();
  dYVector.clear();
  dZVector.clear();
  timeVector.clear();
  
  // the messenger
  primaryMessenger = new CalorPrimaryGeneratorMessenger(this);
  G4cout << " OK!" << G4endl;
}

CalorPrimaryGeneratorAction::~CalorPrimaryGeneratorAction()
{
  delete particleGun;
  
  /*
  delete primaryMessenger;
  if(analysisReaderInstanceFlag)
  {
    delete G4AnalysisReader::Instance();
  }
  if(AngDist) delete AngDist;
  */
}

// --------------------------- set some parameters depending on the primary type
void CalorPrimaryGeneratorAction::SetPrimaryType(G4int type)
{
  typeOfPrimary = type;
  if(typeOfPrimary == 1)   // Particle
  {
    StrTypeOfPrimary = "Particle";
    particleTable = G4ParticleTable::GetParticleTable();
    // gamma is the default particle. Can be changed in the macro by /gun/particle 
    particleName = particleGun->GetParticleDefinition()->GetParticleName();
    primaryParticleDef = particleTable->FindParticle(particleName); 
    particleGun->SetParticleEnergy(particleEnergy);
    PDGencoding = primaryParticleDef->GetPDGEncoding();
    particleTime = 0.;
    
    if(angDistFlag)
      AngDist->ReadFromFile(WkqFile);
  }
  
  else if(typeOfPrimary == 2)  // Isotope
  {
    StrTypeOfPrimary = "Isotope";
    particleGun->SetParticleDefinition(primaryParticleDef);
    particleEnergy   = 0.0*keV;                          // radioisotope at rest    
    particleGun->SetParticleEnergy(particleEnergy);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.0,0.0,1.0)); // for radioisotope
    particleTime = 0.;
  }
  
  else if(typeOfPrimary == 3)  // Phase Space File
  {        
    StrTypeOfPrimary = "Phase Space File";
    // Create (or get) analysis reader
    //if(analysisReaderInstanceFlag)
    //{
      //delete G4AnalysisReader::Instance();
      //analysisReader=0;
    //}
    analysisReader = G4AnalysisReader::Instance();
    analysisReaderInstanceFlag = 1;
    analysisReader->SetVerboseLevel(1);
    analysisReader->SetFileName(phaseSpaceFile);
        
    // ----  GetNtuple()
    // The next line works when the Ntuple is created whithout directory name inside the .root file
    // For the new version of PhaseSpace, the fila is writen as a tree 
    ntupleID = analysisReader->GetNtuple("tree");
        
    if (ntupleID >= 0) 
    {
      G4cout << "tree found at file " << phaseSpaceFile << G4endl;
      analysisReader->SetNtupleIColumn("Reaction", reactionNumber);
      analysisReader->SetNtupleIColumn("Nhits",    Nhits);
      analysisReader->SetNtupleIColumn("PDGchit",  PDGcodeVector);
      analysisReader->SetNtupleIColumn("Volhit",   volumeNbVector); // Targ/Back:1/2 ; Sphere:3                        
      analysisReader->SetNtupleDColumn("Enhit",    enVector);
      analysisReader->SetNtupleDColumn("dXhit",    dXVector);
      analysisReader->SetNtupleDColumn("dYhit",    dYVector);
      analysisReader->SetNtupleDColumn("dZhit",    dZVector);
      analysisReader->SetNtupleDColumn("Thit",     timeVector);
    }
    else 
    {
      G4cout << "Could not find tree at file " << phaseSpaceFile << G4endl;
      G4cout << "Setting type of primary as Particle" << G4endl;
      typeOfPrimary = 1;
      SetPrimaryType(typeOfPrimary);
    }
    // Clear the vectors for accumulation of Ntuple rows
    PDGcodeVector.assign(PDGcodeVector.size(),0);
    PDGcodeVector.clear();
    volumeNbVector.assign(volumeNbVector.size(),0);
    volumeNbVector.clear();
    enVector.assign(enVector.size(),0);
    enVector.clear();
    dXVector.assign(dXVector.size(),0);
    dXVector.clear();
    dYVector.assign(dYVector.size(),0);
    dYVector.clear();
    dZVector.assign(dZVector.size(),0);
    dZVector.clear();
    timeVector.assign(timeVector.size(),0);
    timeVector.clear();
    lastReactionNumber = reactionNumber = 1;
  }
  
  else if(typeOfPrimary == 4)  // Detector Intrinsic Radioactivity
  {
    StrTypeOfPrimary = "Intrinsic Radioactivity";
    G4cout << "Simulation of Instrinsic Radioactivity ..." << G4endl;
    
    G4Material* matScint = detCon->GetDetectorElement()->getMaterialScintillator();
    G4cout << "Material: " << matScint->GetName() << G4endl;
    
    if(matScint->GetName() == "LYSO")
    {  Z_number = 71; A_number = 176; }  // Lu-176 in LSO ~307 Bq/cm3
    else if(matScint->GetName() == "LaBr3")
    {  Z_number = 57; A_number = 138; }  // La-138 in LaBr3
    else
    {
      G4cout << " This scintillator has no radioactive isotope." << G4endl;
      G4cout << " Particle primary type will be used." << G4endl; 
      typeOfPrimary = 1;
      SetPrimaryType(typeOfPrimary);
      return;
    }
    
        
    // Set parameters according the detector shape: Box or Cylinder 
    if(detCon->GetDetectorShape() == "Box") detShape = 1;
    else detShape = 2;  // default shape Cylinder
                  
    if(detShape == 1)
      scintDimX = scintDimY = detCon->GetScintillatorSection();
    else
      scintRadius = detCon->GetScintillatorSection()/2.0;   // diameter/2

    scintDimZ = detCon->GetScintillatorHeight();  // the height in Z is valid for both shapes
    scintPosition = detCon->GetScintillatorPosition();

    particleDirection  = G4ThreeVector(0.,0.,1.);
    // direction has no effect here; emissions will be done from radionuclide decay
    particleGun->SetParticleMomentumDirection(particleDirection);
    // radionuclide at rest inside the scintillator
    particleGun->SetParticleEnergy(0.0);
    G4cout << " OK!" << G4endl;
  }
}


// --------------------------------------------------- the primary generator
void CalorPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //G4cout << G4endl << " PrimaryGeneratorAction::GeneraetePrimaries " << G4endl;
  if(typeOfPrimary == 1)  // Particle
  {
    G4double theta, phi;
    if(!angDistFlag)
    {
      // isotropic distribution for ranges of spherical angles
      // azimuthal angle
      phi   = G4UniformRand()*(phiMax-phiMin)+phiMin;
      // polar angle
      theta = acos(cos(thetaMin)-G4UniformRand()*(cos(thetaMin)-cos(thetaMax)));             
    }
    else  // use sampling for angular distribution
    {
      AngDist->SampleThetaPhi(&theta,&phi);
    }
    
    G4double dir_x = sin(theta)*cos(phi);
    G4double dir_y = sin(theta)*sin(phi);
    G4double dir_z = cos(theta);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(dir_x,dir_y,dir_z));
    PDGencoding = particleGun->GetParticleDefinition()->GetPDGEncoding();
    particleGun->GeneratePrimaryVertex(anEvent);
  }
    
  else if(typeOfPrimary == 2)  // Radiocative Isotope
  {
    primaryParticleDef = G4IonTable::GetIonTable()->GetIon(Z_number,A_number,excitationEnergy);
    particleGun->SetParticleDefinition(primaryParticleDef);
    PDGencoding = particleGun->GetParticleDefinition()->GetPDGEncoding();
    particleGun->GeneratePrimaryVertex(anEvent);
  }
    
  else if(typeOfPrimary == 3)  // Phase Space File
  {
    // Try to read a line of the tree
    if(analysisReader->GetNtupleRow() ) // React, Nhits, PDG, Vol, en, dX, dY, dZ, time
    {
      //G4cout << " --- Primary Generator reading line from input Root file ---" << G4endl; // TEST
      // make a loop of Nhits to generate the primaries
      for(G4int hit = 0; hit < Nhits; hit++)
      {
        //G4cout << "PDG for hit " << hit << ": " << PDGcodeVector[hit] << G4endl;  // TEST
        if(PDGcodeVector[hit] < 1000000000)   // not heavy ion, most frequent
        {
          particleTable = G4ParticleTable::GetParticleTable();
          primaryParticleDef = particleTable->FindParticle(PDGcodeVector[hit]);
        }  
        
        else    // particle is a heavy ion, must be taken from the IonTable
        {
          // Some Ions in excited state can not be used in the particleGun: ERROR
          // The last digit of the PDG for these Ions are different from zero
          // Put a zero in the last digit of these Ions to get them in the ground state            
          if( (PDGcodeVector[hit] - (ceil(PDGcodeVector[hit]/10))*10) > 0)  // test the last digit of PDG
          {
            G4int newPDGencoding = (ceil(PDGcodeVector[hit]/10))*10; // force last digit = 0 (ground state)
            primaryParticleDef = G4IonTable::GetIonTable()->GetIon(newPDGencoding);
          }
          else  // if the last PDG digit = 0, use the original PDG
          {
            primaryParticleDef = G4IonTable::GetIonTable()->GetIon(PDGcodeVector[hit]);
          }
        }   // particle is an Ion
        particleDirection = G4ThreeVector(dXVector[hit],dYVector[hit],dZVector[hit]);
        
        if(volumeNbVector[hit] == 3) // volumeNumber: Target 1; Backing: 2 ; Phase Space Sphere: 3. 
        {
          // position calculation
          particlePosition = radiusOfThePhaseSpaceSphere*particleDirection;
          G4PrimaryVertex* thisVertex = new G4PrimaryVertex(particlePosition,timeVector[hit]);
            
          G4PrimaryParticle* thisPrimaryParticle = new G4PrimaryParticle(primaryParticleDef);
                    
          thisPrimaryParticle->SetMomentumDirection(particleDirection);
          thisPrimaryParticle->SetKineticEnergy(enVector[hit]);
          thisVertex->SetPrimary(thisPrimaryParticle);
          
          anEvent->AddPrimaryVertex(thisVertex); // OLD VERSION          
        }  // if Volume == 3
      }    // for hit (loop for one reaction)
    }      // if GetNtupleRow()
      
    else  // End of the Phase Space file reached. Terminates the present Run.
    {
      endOfPhaseSpaceFileFlag = 1;
      G4cout << G4endl << " End of Phase Space file reached" << G4endl;
      G4cout           << "   This Run will be terminated." << G4endl;
      G4RunManager *RM = G4RunManager::GetRunManager();
      G4int eventsProcessed = RM->GetCurrentRun()->GetNumberOfEvent();
      G4int eventsOrdered   = RM->GetNumberOfEventsToBeProcessed();
      G4cout << "   Number of processed events: " << eventsProcessed << G4endl;
      G4cout << "   Number of ordered events:   " << eventsOrdered << G4endl << G4endl;
      RM->SetNumberOfEventsToBeProcessed(0);
      RM->AbortRun();
      return;
    }
  }
  
  else if(typeOfPrimary == 4)  // Intrinsic Radiation
  {
    // Due to the possibility of changing the detectors parameters with macro commands,
    // this type of primary is dependent on the DetectorConstruction
    // Get the detector's material
    G4Material* matScint = detCon->GetDetectorElement()->getMaterialScintillator();
    if(matScint->GetName() == "LYSO")
    {  Z_number = 71; A_number = 176;}  // Lu-176 in LSO ~307 Bq/cm3
    else if(matScint->GetName() == "LaBr3")
    {  Z_number = 57; A_number = 138;}  // La-138 in LaBr3
    else
    {
      G4cout << "\n This scintillator has no radioactive isotope." << G4endl;
      G4cout << "   This Run will be terminated.\n" << G4endl;
      G4RunManager *RM = G4RunManager::GetRunManager();
      RM->SetNumberOfEventsToBeProcessed(0);
      RM->AbortRun();            
      typeOfPrimary = 1;
      SetPrimaryType(typeOfPrimary);
      return;
    }
    
    // sample one detector
    G4int det = (floor(G4UniformRand()*(detCon->GetNumberOfDetectors())));
    
    // sample a position inside the scintillator
    G4double posX, posY, posZ;
    scintDimZ = detCon->GetScintillatorHeight();  // the height in Z is valid for both shapes
    // The scintillator can be shifted in the Z direction. We take the Z coordinate
    scintPosition = detCon->GetScintillatorPosition();
    G4double posZ_0 = scintPosition.getZ();
    if(detCon->GetDetectorShape() == "Box")
    {
      scintDimX = scintDimY = detCon->GetDetectorSection();          
      posX = -scintDimX/2.+G4UniformRand()*scintDimX;
      posY = -scintDimY/2.+G4UniformRand()*scintDimY;
      posZ = -scintDimZ/2.+G4UniformRand()*scintDimZ + posZ_0;
    }
    else
    {
      scintRadius = detCon->GetDetectorSection()/2.0;
      G4double randomTheta  = 2.0*M_PI*G4UniformRand();          // theta sampling
      G4double randomRadius = scintRadius*sqrt(G4UniformRand()); // radius sampling
      posX = randomRadius*cos(randomTheta);
      posY = randomRadius*sin(randomTheta);
      posZ = -scintDimZ/2.+G4UniformRand()*scintDimZ + posZ_0;   // height sampling
    }
    G4ThreeVector posInside = G4ThreeVector(posX,posY,posZ);
//    G4cout << " Random position:      " << posInside/mm << " mm" << G4endl;
    // rotate and translate to place it at the detector position
//    G4cout << " Detector number:      " << det << G4endl;
//    G4cout << " Detector Angle:       " << detCon->GetDetectorAngle()[det]/deg << " deg" << G4endl;
    G4ThreeVector posNew = posInside.rotateZ(detCon->GetDetectorAngle()[det].getZ());
//    G4cout << " After rotation Psi:   " << posNew/mm << " mm" << G4endl;
    posNew = posNew.rotateY(detCon->GetDetectorAngle()[det].getY());
//    G4cout << " After rotation Theta: " << posNew/mm << " mm" << G4endl;
    posNew = posNew.rotateZ(detCon->GetDetectorAngle()[det].getX());
//    G4cout << " After rotation Phi:   " << posNew/mm << " mm" << G4endl;
    posNew += detCon->GetDetectorPosition()[det];
//    G4cout << " Detector position:    " << detCon->GetDetectorPosition()[det] << G4endl;
//    G4cout << " After translation:    " << posNew/mm << " mm" << G4endl;
    particleGun->SetParticlePosition(posNew);
    primaryParticleDef = G4IonTable::GetIonTable()->GetIon(Z_number,A_number,0.0);
    particleGun->SetParticleDefinition(primaryParticleDef);
    PDGencoding = particleGun->GetParticleDefinition()->GetPDGEncoding();
    particleGun->GeneratePrimaryVertex(anEvent);
  }
}

// ----------------------------------------------------- 
void CalorPrimaryGeneratorAction::SetAngDistDirectory(G4String dir)
{
  WkqDirectory = dir;
  WkqFile      = WkqDirectory+WkqName;
//  SetPrimaryType(1);
}

// ----------------------------------------------------- 
void CalorPrimaryGeneratorAction::SetAngDistName(G4String name)
{
  WkqName      = name;
  WkqFile      = WkqDirectory+WkqName;
  SetPrimaryType(1);
}

// ----------------------------------------------------- 
void CalorPrimaryGeneratorAction::SetAngDistFlag(G4bool flag)
{
  angDistFlag = flag;
  SetPrimaryType(1);
}

// ----------------------------------------------------- 
void CalorPrimaryGeneratorAction::SetPhaseSpaceFileDirectory(G4String dir)
{
  phaseSpaceDirectory = dir;
  phaseSpaceFile      = phaseSpaceDirectory+phaseSpaceName;
  SetPrimaryType(3);
}

// ----------------------------------------------------- 
void CalorPrimaryGeneratorAction::SetPhaseSpaceName(G4String name)
{
  phaseSpaceName      = name;
  phaseSpaceFile      = phaseSpaceDirectory+phaseSpaceName;
  SetPrimaryType(3);
}

// ------------------------------------------------ choose the radioisotope
void CalorPrimaryGeneratorAction::SetIon(G4int Z,G4int A,G4double exEn = 0.)
{
  typeOfPrimary    = 2;                                // ion
  Z_number         = Z;                                // atomic number
  A_number         = A;                                // mass number
  excitationEnergy = exEn;                             // excited state
  particleEnergy = 0.;
}


// ----------------------------------------------------------- print info
void CalorPrimaryGeneratorAction::printPrimaryGeneratorInfo(std::ofstream *output)
{
  *output << G4endl << " --- PrimaryGenerator info for this Run" << G4endl;
  *output     << " Type of primaries:  " << StrTypeOfPrimary << G4endl;
  if(typeOfPrimary == 1)
  {
    particleEnergy  = particleGun->GetParticleEnergy();
    primaryParticleDef = particleGun->GetParticleDefinition();
    *output << " Particle:           " << primaryParticleDef->GetParticleName() << G4endl;
    *output << " Energy:             " << particleEnergy/MeV << " MeV" << G4endl;
    *output << " Source position:    " << particlePosition/mm << " mm" << G4endl;
    if(angDistFlag)
    {
      *output << "\n Particle emitted with PDF for angular distribution\n";
      *output << " File with coefficients:" << G4endl;
      *output << WkqFile << G4endl;
      *output << " Coefficients of the spherical harmonics:" << G4endl;
      for(G4int k = 0; k<=kMax; k++)
        for(G4int q = -k; q<=k; q++)
          *output << " Wkq(" << k << "," << q <<") = " << AngDist->get(k,q) << G4endl;
    }
  }
  
  if(typeOfPrimary == 2)
  {
    primaryParticleDef = G4IonTable::GetIonTable()->GetIon(Z_number,A_number,excitationEnergy);
    particleEnergy  = particleGun->GetParticleEnergy();        
    *output << " Isotope:            " << primaryParticleDef->GetParticleName() << G4endl;
    *output << " Isotope excitation: " << excitationEnergy/MeV << " MeV" << G4endl;
    *output << " Kinetic Energy:     " << particleEnergy/MeV << " MeV" << G4endl;
    *output << " Source position:    " << particlePosition/mm << " mm" << G4endl;
  }
  
  if(typeOfPrimary == 3)
  {
    *output << " File Directory:     " << G4endl << " " << phaseSpaceDirectory << G4endl;
    *output << " File Name:          " << phaseSpaceName  << G4endl;
    if(endOfPhaseSpaceFileFlag)
    {
      *output << " End of Phase Space File reached." << G4endl;
      *output << " All events of the Phase Space File were used." << G4endl;
    }
  }
  
  if(typeOfPrimary == 4)
  {
    *output << " Intrinsic Radiation" << G4endl;
    *output << " Detector material: " << detCon->GetDetectorElement()->getMaterialScintillator()->GetName() << G4endl;
    *output << " Detector shape:    " << detCon->GetDetectorShape() << G4endl;
    if(detCon->GetDetectorShape() == "Box")
      *output << " Detector Side:     " << detCon->GetScintillatorSection() << G4endl;
    else
      *output << " Detector Radius:   " << detCon->GetDetectorSection()/2.0 << G4endl;
    
    *output << " Detector Height:   " << detCon->GetScintillatorHeight() << G4endl;
  }
  *output << G4endl;
  return;
}

//-------------------------------- 
