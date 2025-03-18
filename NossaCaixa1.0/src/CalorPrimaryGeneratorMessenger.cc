// CalorPrimaryGeneratorMessenger

#include "CalorPrimaryGeneratorMessenger.hh"
#include "CalorPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
//#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3Vector.hh"

/////////////////////////////////////////////////////////////////////////////
CalorPrimaryGeneratorMessenger::CalorPrimaryGeneratorMessenger(CalorPrimaryGeneratorAction* PG)
:PrimGen(PG)
{
    // Change Phantom size
    primaryDir = new G4UIdirectory("/primary/");
    primaryDir -> SetGuidance("Commands to change some primary generator parameters");

    // Commands to change angle limits when using particle emission 
    changePhiMinCmd = new G4UIcmdWithADoubleAndUnit("/primary/PhiMin",this);
    changePhiMinCmd->SetDefaultUnit("deg");;
    changePhiMinCmd->SetGuidance("Enter new MINIMUM value for PHI and unit");
    changePhiMinCmd->AvailableForStates(G4State_Idle);

    changePhiMaxCmd = new G4UIcmdWithADoubleAndUnit("/primary/PhiMax",this);
    changePhiMaxCmd->SetDefaultUnit("deg");;
    changePhiMaxCmd->SetGuidance("Enter new MAXIMUM value for PHI and unit");
    changePhiMaxCmd->AvailableForStates(G4State_Idle);

    changeThetaMinCmd = new G4UIcmdWithADoubleAndUnit("/primary/ThetaMin",this);
    changeThetaMinCmd->SetDefaultUnit("deg");;
    changeThetaMinCmd->SetGuidance("Enter new MINIMUM value for THETA and unit");
    changeThetaMinCmd->AvailableForStates(G4State_Idle);

    changeThetaMaxCmd = new G4UIcmdWithADoubleAndUnit("/primary/ThetaMax",this);
    changeThetaMaxCmd->SetDefaultUnit("deg");
    changeThetaMaxCmd->SetGuidance("Enter new MAXIMUM value for THETA and unit");
    changeThetaMaxCmd->AvailableForStates(G4State_Idle);
    
    // Commands to set the Angular Distribution options
    changeAngDistCmd = new G4UIcmdWithABool("/primary/AngularDistribution",this);
    changeAngDistCmd->SetGuidance("Set angular distribution ON: 1 or OFF: 0");
    changeAngDistCmd->SetDefaultValue("0");
    changeAngDistCmd->AvailableForStates(G4State_Idle);
    
    changeKmaxCmd = new G4UIcmdWithAnInteger("/primary/kmax",this);
    changeKmaxCmd->SetGuidance("Change the kmax for coefficients of angular distribution");
    changeKmaxCmd->SetGuidance("Enter new Kmax value for the spherical harmonics coefficients");
    changeKmaxCmd->AvailableForStates(G4State_Idle);
    
    changeAngDistDirectoryCmd = new G4UIcmdWithAString("/primary/AngularDistDirectory", this);
    changeAngDistDirectoryCmd->SetGuidance("Change the DIRECTORY of the file of the Wkq coefficients for Angular Distribution ");
    changeAngDistDirectoryCmd-> SetDefaultValue("./");
    changeAngDistDirectoryCmd-> AvailableForStates(G4State_Idle);
    
    changeAngDistNameCmd = new G4UIcmdWithAString("/primary/AngularDistName", this);
    changeAngDistNameCmd->SetGuidance("Change the NAME of the file of the Wkq coefficients for Angular Distribution");
    changeAngDistNameCmd-> SetDefaultValue("Wkq_examplo.dat");
    changeAngDistNameCmd-> AvailableForStates(G4State_Idle);
    
    // Commands to change some parameters when Phase Space File is used 
    changeRadiusOfSphereCmd = new G4UIcmdWithADoubleAndUnit("/primary/RadiusOfSphere",this);
    changeRadiusOfSphereCmd->SetGuidance("Change the radius of the sphere for Phase Space generation");
    changeRadiusOfSphereCmd->SetDefaultUnit("mm");;
    changeRadiusOfSphereCmd->SetGuidance("Enter new radius value and unit");
    changeRadiusOfSphereCmd->AvailableForStates(G4State_Idle);
    
    
    
    changePhaseSpaceDirectoryCmd = new G4UIcmdWithAString("/primary/PhaseSpaceDirectory", this);
    changePhaseSpaceDirectoryCmd->SetGuidance("Change the DIRECTORY of the Phase Space File");
    changePhaseSpaceDirectoryCmd-> SetDefaultValue("./");
    changePhaseSpaceDirectoryCmd-> AvailableForStates(G4State_Idle);
    
    changePhaseSpaceNameCmd = new G4UIcmdWithAString("/primary/PhaseSpaceName", this);
    changePhaseSpaceNameCmd->SetGuidance("Change the NAME of the Phase Space File");
    changePhaseSpaceNameCmd-> SetDefaultValue("PhSpace_01.root");
    changePhaseSpaceNameCmd-> AvailableForStates(G4State_Idle);
 
    // Command to change the type of primaries
    changePrimaryTypeCmd = new G4UIcmdWithAString("/primary/Type", this);
    changePrimaryTypeCmd->SetGuidance("Change the TYPE of primary: Particle, Isotope, PhaseSpaceFile or IntrinsicRad");
    changePrimaryTypeCmd-> SetDefaultValue("Particle");
    changePrimaryTypeCmd-> AvailableForStates(G4State_Idle);

    changeIonCmd = new G4UIcmdWith3Vector("/primary/Ion", this);
    changeIonCmd->SetGuidance("Change the Ion Z A E");
    changeIonCmd-> AvailableForStates(G4State_Idle);
    
}

/////////////////////////////////////////////////////////////////////////////
CalorPrimaryGeneratorMessenger::~CalorPrimaryGeneratorMessenger()
{
    delete primaryDir; 
    delete changePhaseSpaceDirectoryCmd; 
    delete changePhaseSpaceNameCmd; 
    delete changePrimaryTypeCmd;
    
    delete changeAngDistCmd;
    delete changeKmaxCmd;
    delete changeAngDistDirectoryCmd;
    delete changeAngDistNameCmd;
    
    delete changeRadiusOfSphereCmd;
    delete changeIonCmd;
}

/////////////////////////////////////////////////////////////////////////////
void CalorPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
    if( command == changePhiMinCmd)
    {
        PrimGen->SetPhiMin(changePhiMinCmd->GetNewDoubleValue(newValue));
    }
    else if( command == changePhiMaxCmd)
    {
        PrimGen->SetPhiMax(changePhiMaxCmd->GetNewDoubleValue(newValue));
    }
    else if( command == changeThetaMinCmd)
    {
        PrimGen->SetThetaMin(changeThetaMinCmd->GetNewDoubleValue(newValue));
    }
    else if( command == changeThetaMaxCmd)
    {
        PrimGen->SetThetaMax(changeThetaMaxCmd->GetNewDoubleValue(newValue));
    }
    
    else if( command == changeAngDistCmd)
    {
        PrimGen->SetAngDistFlag(changeAngDistCmd->GetNewBoolValue(newValue));
    }
    else if( command == changeKmaxCmd)
    {
        PrimGen->SetKmax(changeKmaxCmd->GetNewIntValue(newValue));
    }
    else if (command == changeAngDistDirectoryCmd)
    {
        PrimGen->SetAngDistDirectory(newValue);
    }
    else if (command == changeAngDistNameCmd)
    {
        PrimGen->SetAngDistName(newValue);
    }

    
    else if( command == changeRadiusOfSphereCmd)
    {
        PrimGen->SetRadiusOfSphere(changeRadiusOfSphereCmd->GetNewDoubleValue(newValue));
    }
    
    else if (command == changePhaseSpaceDirectoryCmd)
    {
        PrimGen->SetPhaseSpaceFileDirectory(newValue);
    }
    else if (command == changePhaseSpaceNameCmd)
    {
        PrimGen->SetPhaseSpaceName(newValue);
    }
    else if (command == changePrimaryTypeCmd)
    {
        G4int primaryType = 1;
        if(newValue == "Isotope")
            primaryType = 2;
        else if(newValue == "PhaseSpaceFile")
            primaryType = 3;
        else if(newValue == "IntrinsicRad")
            primaryType = 4;
        PrimGen->SetPrimaryType(primaryType);
    }
    else if (command == changeIonCmd)
    {
        G4ThreeVector ZAE = changeIonCmd->GetNew3VectorValue(newValue);
        G4int Z = G4int(ZAE.x());
        G4int A = G4int(ZAE.y());
        PrimGen->SetIon(Z,A,ZAE.z());
    }

}
