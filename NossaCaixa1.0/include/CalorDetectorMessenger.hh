//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: CalorDetectorMessenger.hh
// -------------------------------------------------------------------

#ifndef CalorDetectorMessenger_h
#define CalorDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CalorDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
//class G4UIcmdWithAnInteger;
//class G4UIcmdWith3VectorAndUnit;

class CalorDetectorMessenger: public G4UImessenger
{
  public:
    CalorDetectorMessenger(CalorDetectorConstruction* );
   ~CalorDetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String);

 private:
    CalorDetectorConstruction *DetCon;
    
    G4UIdirectory        *detectorsDir;
    G4UIdirectory        *chamberDir;

    G4UIcmdWithAString         *changeFileDirectoryCmd; 
    G4UIcmdWithAString         *changeDescriptorNameCmd;
    G4UIcmdWithAString         *changeDetectorMaterialCmd;
    G4UIcmdWithAString         *changeDetectorShapeCmd;
 //   G4UIcmdWithADoubleAndUnit  *changeDetectorLengthCmd;
 //   G4UIcmdWithADoubleAndUnit  *changeDetectorSectionCmd;
    G4UIcmdWithADoubleAndUnit  *changeRadialDislocationCmd;
    G4UIcmdWithAString         *changeChamberMaterialCmd;    
    G4UIcmdWithADoubleAndUnit  *changeChamberThicknessCmd; 
    G4UIcmdWithADoubleAndUnit  *changeChamberOuterRadiusCmd; 
    G4UIcmdWithAString         *changeTargetHolderMaterialCmd;
};

#endif





