// $Id: CalorDetectorMessenger.cc

#include "CalorDetectorMessenger.hh"
#include "CalorDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
//#include "G4UIcmdWithAnInteger.hh"
//#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4SystemOfUnits.hh"

//------------------------------------------------------

CalorDetectorMessenger::CalorDetectorMessenger(CalorDetectorConstruction * DC):
DetCon(DC)
{

    // Change Detectors parameters
    detectorsDir = new G4UIdirectory("/detectors/");
    detectorsDir -> SetGuidance("Commands to change the detectors parameters");
    
    changeFileDirectoryCmd = new G4UIcmdWithAString("/detectors/FileDirectory", this);
    changeFileDirectoryCmd->SetGuidance("Change the DIRECTORY of the File");
    changeFileDirectoryCmd->SetDefaultValue("./");
    changeFileDirectoryCmd->AvailableForStates(G4State_Idle);
    
    changeDescriptorNameCmd = new G4UIcmdWithAString("/detectors/FileName", this);
    changeDescriptorNameCmd->SetGuidance("Change the NAME of the File");
    changeDescriptorNameCmd->SetDefaultValue("lplace-radial.dat");
    changeDescriptorNameCmd->AvailableForStates(G4State_Idle);
    
    changeDetectorMaterialCmd = new G4UIcmdWithAString("/detectors/Material", this);
    changeDetectorMaterialCmd->SetGuidance("Change the Material (LYSO, LaBr3 or GAGG)");
    changeDetectorMaterialCmd->SetDefaultValue("LaBr3");
    changeDetectorMaterialCmd->AvailableForStates(G4State_Idle);
    
    changeDetectorShapeCmd = new G4UIcmdWithAString("/detectors/Shape", this);
    changeDetectorShapeCmd->SetGuidance("Change the shape (Box or Cylinder)");
    changeDetectorShapeCmd->SetDefaultValue("Box");
    changeDetectorShapeCmd->AvailableForStates(G4State_Idle);

/*    
    changeDetectorLengthCmd = new G4UIcmdWithADoubleAndUnit("/detectors/Length", this);
    changeDetectorLengthCmd->SetGuidance("Change the Detectors Length");
    changeDetectorLengthCmd->SetDefaultValue(50.0*mm);
    changeDetectorLengthCmd->AvailableForStates(G4State_Idle);

    changeDetectorSectionCmd = new G4UIcmdWithADoubleAndUnit("/detectors/Section", this);
    changeDetectorSectionCmd->SetGuidance("Change the Detectors Side or Diameter");
    changeDetectorSectionCmd->SetDefaultValue(25.0*mm);
    changeDetectorSectionCmd->AvailableForStates(G4State_Idle);
*/    
    changeRadialDislocationCmd = new G4UIcmdWithADoubleAndUnit("/detectors/Dislocation", this);
    changeRadialDislocationCmd->SetGuidance("Change the Radial Dislocation (value and unit)");
    changeRadialDislocationCmd->SetDefaultValue(0.0*mm);
    changeRadialDislocationCmd->AvailableForStates(G4State_Idle);
    
    // Change Chamber parameters
    chamberDir = new G4UIdirectory("/chamber/");
    chamberDir -> SetGuidance("Commands to change the chamber parameters");

    changeChamberMaterialCmd = new G4UIcmdWithAString("/chamber/Material", this);
    changeChamberMaterialCmd->SetGuidance("Change the Material (G4_Al or G4_Galactic (vacuum))");
    changeChamberMaterialCmd->SetDefaultValue("G4_Al");
    changeChamberMaterialCmd->AvailableForStates(G4State_Idle);
    
    changeChamberThicknessCmd = new G4UIcmdWithADoubleAndUnit("/chamber/Thickness", this);
    changeChamberThicknessCmd->SetGuidance("Change the Chamber Thickness");
    changeChamberThicknessCmd->SetDefaultValue(10.0*mm);
    changeChamberThicknessCmd->AvailableForStates(G4State_Idle);
    
    changeChamberOuterRadiusCmd = new G4UIcmdWithADoubleAndUnit("/chamber/Radius", this);
    changeChamberOuterRadiusCmd->SetGuidance("Change the Chamber External Radius");
    changeChamberOuterRadiusCmd->SetDefaultValue(240.0*mm);
    changeChamberOuterRadiusCmd->AvailableForStates(G4State_Idle);
    
    changeTargetHolderMaterialCmd = new G4UIcmdWithAString("/chamber/TargetHolderMaterial", this);
    changeTargetHolderMaterialCmd->SetGuidance("Change the Material (G4_Al, G4_Cu or G4_Galactic)");
    changeTargetHolderMaterialCmd->SetDefaultValue("G4_Cu");
    changeTargetHolderMaterialCmd->AvailableForStates(G4State_Idle);
    
}

//------------------------------------------------------
CalorDetectorMessenger::~CalorDetectorMessenger()
{
    delete changeDescriptorNameCmd;
    delete changeDetectorMaterialCmd;
    delete changeDetectorShapeCmd;
//    delete changeDetectorLengthCmd;
//    delete changeDetectorSectionCmd;
    delete changeRadialDislocationCmd;
    delete changeChamberThicknessCmd;
    delete changeChamberOuterRadiusCmd;
    delete changeTargetHolderMaterialCmd;
    delete changeFileDirectoryCmd; 
    delete chamberDir;    
    delete detectorsDir; 
}


//------------------------------------------------------
void CalorDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
    if (command == changeFileDirectoryCmd)
    {
        DetCon->SetDetectorsFileDirectory(newValue);
    }
    else if (command == changeDescriptorNameCmd)
    {
        DetCon->SetDetectorsDescriptorName(newValue);
    }
    else if (command == changeDetectorMaterialCmd)
    {
        DetCon->SetDetectorsMaterial(newValue);
    }
    else if (command == changeDetectorShapeCmd)
    {
        DetCon->SetDetectorsShape(newValue);
    }

//    else if( command == changeDetectorLengthCmd)
//    {
//        DetCon->SetDetectorsLength(changeDetectorLengthCmd->GetNewDoubleValue(newValue));
//    }
//    else if( command == changeDetectorSectionCmd)
//    {
//        DetCon->SetDetectorsSection(changeDetectorSectionCmd->GetNewDoubleValue(newValue));
//    }

else if( command == changeRadialDislocationCmd)
    {
        DetCon->SetRadialDislocation(changeRadialDislocationCmd->GetNewDoubleValue(newValue));
    }
    else if (command == changeChamberMaterialCmd)
    {
        DetCon->SetChamberMaterial(newValue);
    }
    else if( command == changeChamberThicknessCmd)
    {
        DetCon->SetChamberThickness(changeChamberThicknessCmd->GetNewDoubleValue(newValue));
    }
    else if( command == changeChamberOuterRadiusCmd)
    {
        DetCon->SetChamberOuterRadius(changeChamberOuterRadiusCmd->GetNewDoubleValue(newValue));
    }
    else if (command == changeTargetHolderMaterialCmd)
    {
        DetCon->SetTargetHolderMaterial(newValue);
    }
    
}


