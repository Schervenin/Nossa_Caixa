
// $Id: Detector.hh
// GEANT4 tag $Name:  $

// Version 3.0: option of cylindrical geometry

// Detector: scintillator with square section. and a cover
// to protect the scintillator and/or to reflect the light


#ifndef DetectorElement_H
#define DetectorElement_H 1

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4VisAttributes;
class G4Box;
class G4Tubs;

class ofstream;

class DetectorElement
{
  public:
    DetectorElement();
    DetectorElement(G4String material, G4String shape, G4String nameDef);
   ~DetectorElement();
     
  private:
    void DefineMaterials();  // define the materials
    void DefineColors();
    void BuildDetector();
    
  private:
    G4String name;           // a name for this detector
    G4String scintiName;     // a name for the scintillator volume
    G4String type;           // Box or Cylinder
    
// Note on the geometry convention:
// There is simetry around the Z axis
// "Section" and "Thickness" are used for the X-Y plane.  
// "Height"                   is used for the Z direction
// "Section" is used for side (for Box type) or diameter (for Cylinder type)

    // --------------------------- Detector
    // --- A Mother Volume to contain the whole detector element
    G4Box*           detMotherBox_solid = nullptr;
    G4Tubs*          detMotherCyl_solid = nullptr;
    G4LogicalVolume* detMother_log = nullptr;
    G4double         detMotherSection, detMotherHeight;

/* / --- A cap to cover the detector element (aluminum)
    G4Material*      MaterialCover = nullptr;
// ------ The side
    G4Box*           coverSideBox_solid = nullptr;
    G4Tubs*          coverSideCyl_solid = nullptr;
    G4LogicalVolume* coverSide_log = nullptr;
    G4double         coverSideSection, coverSideHeight;
    G4double         coverSideThickness;
    G4ThreeVector    coverSidePosition;
// ------ The front
    G4Box*           coverFrontBox_solid = nullptr;
    G4Tubs*          coverFrontCyl_solid = nullptr;
    G4LogicalVolume* coverFront_log = nullptr;
    G4double         coverFrontSection;
    G4double         coverFrontHeight;
    G4ThreeVector    coverFrontPosition;
*/

// --- The Reflector (BaSO4) only at the sides
    G4Material*      MaterialReflector = nullptr;    
    G4Box*           reflectorSideBox_solid = nullptr;
    G4Tubs*          reflectorSideCyl_solid = nullptr;
    G4LogicalVolume* reflectorSide_log = nullptr;
    G4double         reflectorSideSection, reflectorSideHeight;
    G4double         reflectorSideThickness;
    G4ThreeVector    reflectorSidePosition;
    
// --- The Teflon
    G4Material*      MaterialTeflon;
// ------ The side
    G4Box*           teflonSideBox_solid = nullptr;
    G4Tubs*          teflonSideCyl_solid = nullptr;
    G4LogicalVolume* teflonSide_log = nullptr;
    G4double         teflonSideSection, teflonSideHeight;
    G4double         teflonSideThickness;
    G4ThreeVector    teflonSidePosition;
// ------ The front
    G4Box*           teflonFrontBox_solid = nullptr;
    G4Tubs*          teflonFrontCyl_solid = nullptr;
    G4LogicalVolume* teflonFront_log = nullptr;
    G4double         teflonFrontSection;
    G4double         teflonFrontHeight;
    G4ThreeVector    teflonFrontPosition;

// --- The Scintillator
    G4Material*       MaterialScintillator = nullptr;
    G4Box*            boxScintillator_solid = nullptr;
    G4Tubs*           cylScintillator_solid = nullptr;
    G4LogicalVolume*  scintillator_log = nullptr;
    G4double          scintillatorSection,scintillatorHeight;
    G4ThreeVector     scintillatorPosition;

    
  public:
    G4String          getName()                  { return name;}
    G4String          getScintillatorName()      { return scintiName;}
    
    G4LogicalVolume*  getDetectorPointer()       { return detMother_log; }    
    G4double          getDetectorSection()       { return detMotherSection;}
    G4double          getDetectorHeight()        { return detMotherHeight;}

    G4LogicalVolume*  getScintiPointer()         { return scintillator_log; }
    
    G4Material*       getMaterialScintillator()  { return MaterialScintillator;}
    G4double          getScintillatorSection()   { return scintillatorSection;}    
    G4double          getScintillatorHeight()    { return scintillatorHeight;}        
    G4ThreeVector     getScintillatorPosition()  { return scintillatorPosition;}

//    void setScintillatorMaterial(G4String mat);
//    void setScintillatorShape(G4String shape);
//    void setScintillatorHeight(G4double height);
//    void setScintillatorSection(G4double section);
    
    
    void printInfo(std::ofstream *outputInfo);
        
   private:
// Colors
   G4VisAttributes* blue_scintillator = nullptr;
   G4VisAttributes* yellow_reflector = nullptr;
   G4VisAttributes* white_teflon = nullptr;
   //G4VisAttributes* red_cover = nullptr;
};

#endif
