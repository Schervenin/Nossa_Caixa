
// DetectorElement.cc
// Class that builds a gamma ray detector set

// This class has a method to produce the energy spectrum with resolution
// given by a rough calibration depending on the material 

// Geometry:  parallelepiped with square section
//            cover made of aluminun with different side and frontal thicknesses
// Materials: LaBr3, LYSO, NaI, HPGe


// Author: Maurício Moralles

// version 3.1 - more details: reflector, teflon, cover

// version 3.0 - date 23/06/2019
//              Inclusion of the option for cylindric shape

// version 1.0 - date:   31/10/2017


#include "DetectorElement.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"   

#include <fstream>
using namespace std;


DetectorElement::DetectorElement(G4String mat = "LYSO", G4String shape = "Box", G4String nameDef = "D_")
//:detMotherBox_solid(0),detMotherCyl_solid(0),detMother_log(0),
//MaterialCover(0),
//coverSideBox_solid(0),coverSideCyl_solid(0),coverSide_log(0),
//coverFrontBox_solid(0),coverFrontCyl_solid(0),coverFront_log(0),
//MaterialReflector(0),
//reflectorSideBox_solid(0),reflectorSideCyl_solid(0),reflectorSide_log(0),
//MaterialTeflon(0),
//teflonSideBox_solid(0),teflonSideCyl_solid(0),teflonSide_log(0),
//teflonFrontBox_solid(0),teflonFrontCyl_solid(0),teflonFront_log(0),
//MaterialScintillator(0),
//boxScintillator_solid(0),cylScintillator_solid(0),scintillator_log(0),
//blue_scintillator(0),yellow_reflector(0),white_teflon(0),red_cover(0)
{
  // Initial definitions
  DefineMaterials();
  DefineColors();
  
  MaterialScintillator = G4Material::GetMaterial(mat);
  if(!MaterialScintillator)
  {
    G4cout << "Material " << mat << " not available." << G4endl;
    G4cout << "LaBr3 will be used." << G4endl;
    MaterialScintillator = G4Material::GetMaterial("LaBr3");
  }
  name = nameDef;
  type = shape;

  // ------------------------------------ Elements dimensions
// Note on the geometry convention:
// There is simetry around the Z axis
// "Section" and "Thickness" are used for the X-Y plane.  
// "Height"                   is used for the Z direction 
// "Section" is used for side (for Box type) or diameter (for Cylinder type)
  
  // ------ Scintillator
  // The next line is a name for the scintilator volume to be used to identify the sensitive
  // volume in the SteppingAction
  scintiName = nameDef + "Scinti";  // name to be given to the scintillator volume 
  scintillatorSection  = 12.4*mm;    // side/diameter of the scintillator  FIXED
  scintillatorHeight   = 40.0*mm;    // height of the scintillator         FIXED
  
  // ------ Reflector
  reflectorSideThickness = 0.4*mm;          // reflector (BaSO4)    FIXED
  reflectorSideSection   = scintillatorSection + 2.0*reflectorSideThickness;
  reflectorSideHeight    = scintillatorHeight;
  
  // ------ Teflon side
  teflonSideThickness = 0.2*mm;            // PTFE tape             FIXED
  teflonSideSection   = reflectorSideSection + 2.0*teflonSideThickness;
  teflonSideHeight    = scintillatorHeight;
  // ------ Teflon front
  teflonFrontSection   = teflonSideSection;
  teflonFrontHeight = 0.2*mm;        // teflon front             FIXED

  /* / ------ Cover (aluminum)
  // ------ Cover Side (< = 0 = >) modified to "no cover"
  coverSideThickness  = 0.0*mm;            // cover side (aluminum) FIXED
  coverSideSection    = scintillatorSection;           // cover side (aluminum) FIXED
  coverSideHeight     = scintillatorHeight+teflonFrontHeight; 
  // ------ Cover Front
  coverFrontSection   = coverSideSection;  // cover front
  coverFrontHeight    = 0.0*mm;            // cover front           FIXED */
    
  // Detector Mother Volume - has the size of the complete cover
  //detMotherSection = coverSideSection;
 // detMotherHeight  = coverSideHeight + coverFrontHeight;
 
 // Detector Mother Volume
  detMotherSection = teflonSideSection;
  detMotherHeight  = scintillatorHeight + teflonFrontHeight;
  
  // --- Positions relative to the Mother Volume -------------------------------
  // The X and Y coordinates are always 0.0 for all daughter volumes
  //  --- CoverSide
  //G4double posZ = detMotherHeight/2.0-coverSideHeight/2.0; // value of Z for Cover Side
 // G4double posZ = detMotherHeight/2.0; // value of Z for Cover Side
 // coverSidePosition = G4ThreeVector(0.0,0.0,posZ);  
  // ---CoverFront
  //posZ = -detMotherHeight/2.0+coverFrontHeight/2.0;     // value of Z for Cover Front
  //coverFrontPosition = G4ThreeVector(0.0,0.0,posZ);        
 
  // --- Scintillator
  G4double posZ = detMotherHeight/2.0-scintillatorHeight/2.0;        // value of Z for Scintillator
  scintillatorPosition = G4ThreeVector(0.0,0.0,posZ);

  // --- ReflectorSide
  posZ = detMotherHeight/2.0-reflectorSideHeight/2.0;
  reflectorSidePosition = G4ThreeVector(0.0,0.0,posZ);
 
  // --- TeflonSide
  posZ = detMotherHeight/2.0-teflonSideHeight/2.0;
  teflonSidePosition = G4ThreeVector(0.0,0.0,posZ);
  // --- TeflonFront
  //posZ = -detMotherHeight/2.0+coverFrontHeight+teflonFrontHeight/2.0;
  posZ = -detMotherHeight/2.0+teflonFrontHeight/2.0;
  teflonFrontPosition = G4ThreeVector(0.0,0.0,posZ);

      
  BuildDetector();
  
  // Screen output
  
  G4cout << G4endl;
  G4cout << "--- Detector properties ---" << G4endl;
  G4cout << "--- Detector Type:      " << type << G4endl;  
  G4cout << "--- Scintillator ---" << G4endl;  
  G4cout << " Scintillator material: " << MaterialScintillator->GetName() << G4endl;
 if(type == "Box")
 {
  G4cout << " Section side:          " << scintillatorSection/mm << " mm" << G4endl;
 }
 else
 {
  G4cout << " Section diameter:      " << scintillatorSection/mm << " mm" << G4endl;
 }
  G4cout << " Height:                " << scintillatorHeight/mm << " mm" << G4endl;

  G4cout << "--- Reflector ---" << G4endl;
  G4cout << " Reflector material:    " << MaterialReflector->GetName() << G4endl;
  G4cout << " Thickness:             " << reflectorSideThickness/mm << " mm" << G4endl;

  G4cout << "--- Teflon ---" << G4endl;
  G4cout << " Teflon material:       " << MaterialTeflon->GetName() << G4endl;
  G4cout << " Side Thickness:        " << teflonSideThickness/mm << " mm" << G4endl;
  G4cout << " Front Thickness:       " << teflonFrontHeight/mm << " mm" << G4endl;
  
 /* G4cout << "--- Cover ---" << G4endl;
  G4cout << " Cover material:        " << MaterialCover->GetName() << G4endl;
  G4cout << " Side  thickness:       " << coverSideThickness/mm << " mm" << G4endl;
  G4cout << " Front thickness:       " << coverFrontHeight/mm << " mm" << G4endl;
*/
  G4cout << " --- Complete detector ---" << G4endl;
 if(type == "Box")
 {
  G4cout << " Section side:          " << detMotherSection/mm << " mm" << G4endl;
 }
 else
 {
  G4cout << " Section diameter:      " << detMotherSection/mm << " mm" << G4endl;    
 }
  G4cout << " Height:                " << detMotherHeight/mm << " mm" << G4endl;
  G4cout << G4endl;
  
}

DetectorElement::DetectorElement()
{
  if(detMotherBox_solid)     delete detMotherBox_solid;
  if(detMotherCyl_solid)     delete detMotherCyl_solid;
  
 /* if(MaterialCover)          delete MaterialCover;
  if(coverSideBox_solid)     delete coverSideBox_solid;
  if(coverSideCyl_solid)     delete coverSideCyl_solid;
  if(coverSide_log)          delete coverSide_log;
  if(coverFrontBox_solid)    delete coverFrontBox_solid;
  if(coverFrontCyl_solid)    delete coverFrontCyl_solid;
  if(coverFront_log)         delete coverFront_log;
*/

  if(MaterialReflector)      delete MaterialReflector;
  if(reflectorSideBox_solid) delete reflectorSideBox_solid;
  if(reflectorSideCyl_solid) delete reflectorSideCyl_solid;
  if(reflectorSide_log)      delete reflectorSide_log;

  if(MaterialTeflon)         delete MaterialTeflon;
  if(teflonSideBox_solid)    delete teflonSideBox_solid;
  if(teflonSideCyl_solid)    delete teflonSideCyl_solid;
  if(teflonSide_log)         delete teflonSide_log;
  if(teflonFrontBox_solid)   delete teflonFrontBox_solid;
  if(teflonFrontCyl_solid)   delete teflonFrontCyl_solid;
  if(teflonFront_log)        delete teflonFront_log;
  
  if(MaterialScintillator)   delete MaterialScintillator;  
  if(boxScintillator_solid)  delete boxScintillator_solid;
  if(cylScintillator_solid)  delete cylScintillator_solid;
  if(scintillator_log)       delete scintillator_log;
}

// ---  Destructor
DetectorElement::~DetectorElement(){;}

// --------------------------------------------------- DefineMaterials
void DetectorElement::DefineMaterials()
{
    G4NistManager* NistMan = G4NistManager::Instance();
    G4bool isotopes = false;
  
//    if(!G4Material::GetMaterial("G4_SODIUM_IODIDE",false))   // scintillator
//      NistMan->FindOrBuildMaterial("G4_SODIUM_IODIDE",isotopes);
//    if(!G4Material::GetMaterial("G4_Pb",false))             
//      NistMan->FindOrBuildMaterial("G4_Pb",isotopes);
    //if(!G4Material::GetMaterial("G4_Al",false))                // for cover
      NistMan->FindOrBuildMaterial("G4_Al",isotopes);
    if(!G4Material::GetMaterial("G4_BARIUM_SULFATE",false))    // for reflector
      NistMan->FindOrBuildMaterial("G4_BARIUM_SULFATE",isotopes);    
//    if(!G4Material::GetMaterial("G4_Ge",false))              // HPGe
//      NistMan->FindOrBuildMaterial("G4_Ge",isotopes);
    if(!G4Material::GetMaterial("G4_TEFLON",false))            // Teflon tape
      NistMan->FindOrBuildMaterial("G4_TEFLON",isotopes);
//    if(!G4Material::GetMaterial("G4_POLYCARBONATE",false))   //  option for cover
//      NistMan->FindOrBuildMaterial("G4_POLYCARBONATE",isotopes);
    if(!G4Material::GetMaterial("G4_AIR",false))              // for mother volume
      NistMan->FindOrBuildMaterial("G4_AIR",isotopes);
//    if(!G4Material::GetMaterial("G4_Galactic",false))       // vacuum
//      NistMan->FindOrBuildMaterial("G4_Galactic",isotopes);

// --------------------------- More materials not defined by NIST      
// LaBr3 (scintillator)
// LYSO  (scintillator)
// GAGG  (scintillator)
//----------------------------------------------------------- ELEMENTS
  G4double a;  // atomic mass
  G4double z;  // atomic number
  
  a = 79.904*g/mole;
  z = 35.;
  G4Element* elBr = new G4Element("Bromine","Br",z,a);
  
  a = 138.9055*g/mole;
  z = 57.;
  G4Element* elLa = new G4Element("Lanthanum","La",z,a);

  a = 15.9994*g/mole;
  z = 8.;
  G4Element* elO = new G4Element("Oxygen","O",z,a);
  
  a = 28.0855*g/mole;
  z = 14.;
  G4Element* elSi = new G4Element("Silicon","Si",z,a);

  a = 88.90585*g/mole;
  z = 39.;
  G4Element* elY = new G4Element("Yttrium","Y",z,a);

  a = 174.967*g/mole;
  z = 71.;
  G4Element* elLu = new G4Element("Lutetium","Lu",z,a);
  
  a = 26.9815*g/mole;
  z = 13.;
  G4Element* elAl = new G4Element("Aluminum","Al",z,a);
  
  a = 69.723*g/mole;
  z = 31.;
  G4Element* elGa = new G4Element("Galium","Ga",z,a);  
  
  a = 157.25*g/mole;
  z = 64.;
  G4Element* elGd = new G4Element("Gadolinium","Gd",z,a);  
  
  
//------------------------------------------------------ MATERIALS
  G4double density;  //, massfraction; 
  G4int    nelements, natoms;
  
  // ------------------------------------ material: LaBr3
  // LaBr3: define now if not yet defined.  false: turn off "warning" if already defined
  if(!G4Material::GetMaterial("LaBr3",false))
  {
      density = 5.08*g/cm3;  // from Saint-Gobain tables
      G4Material *LaBr3 = new G4Material("LaBr3", density, nelements = 2);
      LaBr3->AddElement(elLa, natoms = 1);
      LaBr3->AddElement(elBr, natoms = 3);
      G4cout << "Defining material LaBr3." << G4endl;  
  }
  
  // ------------------------------------ material: LYSO -> Lu(1.9)Y(0.1)Si(1)O(5)
  // Composition from PROTEUS LYSO manufacturer
  if(!G4Material::GetMaterial("LYSO",false))
  {  
      density = 7.2*g/cm3;
      G4Material *LYSO = new G4Material("LYSO", density, nelements = 4);
      LYSO->AddElement(elLu,19);
      LYSO->AddElement(elY,1);
      LYSO->AddElement(elSi,10);
      LYSO->AddElement(elO,50);
      G4cout << "Defining material LYSO." << G4endl;    
  }
  // ------------------------------------ material: GAGG -> Ce:Gd3 Al2 Ga3 O12
  if(!G4Material::GetMaterial("GAGG",false))
  {  
      density = 6.63*g/cm3;
      G4Material *GAGG = new G4Material("GAGG", density, nelements = 4);
      GAGG->AddElement(elGd,3);
      GAGG->AddElement(elAl,2);
      GAGG->AddElement(elGa,3);
      GAGG->AddElement(elO,12);
      G4cout << "Defining material GAGG." << G4endl;    
  }  
}

// ----------------------------------------------------------- Define Colors
void DetectorElement::DefineColors()
{
// --- Colors

// G4Colour( R , G , B , t) t=transparency: 1 = opaque, 0 = total transparency
  blue_scintillator =  new G4VisAttributes(G4Colour(0.3,0.4,0.7,1.0));  // 1.0 = no transparency = opaque
  yellow_reflector  =  new G4VisAttributes(G4Colour(1.0,0.8,0.0,0.4));; // 0.4 = transparency
  white_teflon      =  new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.2));; // 0.2 = transparency
 // red_cover        =  new G4VisAttributes(G4Colour(0.8,0.1,0.1,0.2));; // 0.2 = transparency
  
  blue_scintillator->SetVisibility(true);
  blue_scintillator->SetForceSolid(true);  // blue solid

  yellow_reflector->SetVisibility(true);
  yellow_reflector->SetForceSolid(true);  // yellow solid
  
  white_teflon->SetVisibility(true);
  white_teflon->SetForceSolid(true);      // white solid
  
  //red_cover->SetVisibility(true);
  //red_cover->SetForceSolid(true);        // gray solid  
}


// ------------------------------------------------ Build the DetectorElement
void DetectorElement::BuildDetector()
{
  G4cout << "Building the " << name << " detector of type " << type << G4endl << G4endl;
  
  // --- Mother volume for detector element. Same size as the cover ---------------
  G4Material* material = G4Material::GetMaterial("G4_AIR");  // for the mother volume  
  G4cout << "Mother volume" 
         << "       section:  " << detMotherSection 
         << "       height:   " << detMotherHeight << G4endl;
  G4cout << "       material: " << material->GetName() << G4endl;                                       

  
  G4String nameTemp = name + "Det_sol";
  if(type == "Box")
  {
    detMotherBox_solid = new G4Box(nameTemp,detMotherSection/2.0,detMotherSection/2.0,detMotherHeight/2.0);
    nameTemp = name + "Det_log";
    detMother_log = new G4LogicalVolume(detMotherBox_solid,material,nameTemp,0,0,0);    
  }
  else  // not Box, assumed cylinder
  {
    G4double outerRadius  = detMotherSection/2.0;
    G4double innerRadius  = 0.0;
    G4double initialAngle = 0.0*deg;
    G4double finalAngle   = 360.0*deg;
    detMotherCyl_solid = new G4Tubs(nameTemp,innerRadius,outerRadius,detMotherHeight/2.0,initialAngle,finalAngle);
    nameTemp      = name + "Det_log";
    detMother_log = new G4LogicalVolume(detMotherCyl_solid,material,nameTemp,0,0,0);
  }
 

  // --- Scintillator -------------------------------------------------------
  G4cout << "Inserting the Scintillator " << G4endl 
         << "          section:  " << scintillatorSection 
         << "          height:   " << scintillatorHeight << G4endl;
  G4cout << "          material: " << MaterialScintillator->GetName() << G4endl;        

  nameTemp = scintiName + "_sol";
  if(type == "Box")
  {
    boxScintillator_solid = new G4Box(nameTemp,scintillatorSection/2.0,scintillatorSection/2.0,scintillatorHeight/2.0);
    nameTemp = scintiName + "_log";
    scintillator_log = new G4LogicalVolume(boxScintillator_solid,MaterialScintillator,nameTemp,0,0,0);
  }
  else  // assumed cylinder
  {
    G4double outerRadius  = scintillatorSection/2.0;
    G4double innerRadius  = 0.0;
    G4double initialAngle = 0.0*deg;
    G4double finalAngle   = 360.0*deg;
    cylScintillator_solid = new G4Tubs(nameTemp,innerRadius,outerRadius,scintillatorHeight/2.0,initialAngle,finalAngle);
    nameTemp = scintiName + "_log";    
    scintillator_log = new G4LogicalVolume(cylScintillator_solid,MaterialScintillator,nameTemp,0,0,0);    
  }
 
 
  // --- Reflector -----------------------------------------------------------
  MaterialReflector = G4Material::GetMaterial("G4_BARIUM_SULFATE");  
  G4cout << "Inserting the Reflector " << G4endl 
         << "          section : " << reflectorSideSection 
         << "          height:   " << reflectorSideHeight << G4endl;
  G4cout << "          material: " << MaterialReflector->GetName() << G4endl;        


  if(type == "Box")
  {
    // --- Reflector side
    nameTemp = name + "ReflectorS_sol";
    reflectorSideBox_solid = new G4Box(nameTemp,reflectorSideSection/2.0,reflectorSideSection/2.0,reflectorSideHeight/2.0);
    // The reflector side is a box with with a hollow obtained by subtraction of volumes
    
    // Reflector side hollow definition
    G4double reflectorSideHollowSection = reflectorSideSection - 2.0*reflectorSideThickness;
    G4double reflectorSideHollowHeight  = reflectorSideHeight;
    G4Box *reflectorHollow_solid = new G4Box("ReflecHollow_sol",reflectorSideHollowSection/2.0,reflectorSideHollowSection/2.0,reflectorSideHollowHeight/2.0);
    
    // Now the final solid as a subtraction
    G4SubtractionSolid* reflectorFinal_solid = new G4SubtractionSolid("ReflecSide_sol",reflectorSideBox_solid,reflectorHollow_solid);
    reflectorSide_log = new G4LogicalVolume(reflectorFinal_solid,MaterialReflector,"ReflecSide_log",0,0,0);
  }
  else  // assumed cylinder
  {
    // --- Reflector side
    nameTemp = name + "ReflectorS_sol";
    G4double outerRadius  = reflectorSideSection/2.0;
    G4double innerRadius  = 0.0;
    G4double initialAngle = 0.0*deg;
    G4double finalAngle   = 360.0*deg;
    reflectorSideCyl_solid = new G4Tubs(nameTemp,innerRadius,outerRadius,reflectorSideHeight/2.0,initialAngle,finalAngle);

    // Reflector side hollow definition
    outerRadius  = reflectorSideSection/2.0 - reflectorSideThickness;
    G4double reflectorSideHollowHeight = reflectorSideHeight;
    G4Tubs* reflectorHollow_solid = new G4Tubs("ReflecHollow_sol",innerRadius,outerRadius,reflectorSideHollowHeight/2.0,initialAngle,finalAngle);

    // Now the final solid as a subtraction
    G4SubtractionSolid* reflectorFinal_solid = new G4SubtractionSolid("ReflecSide_sol",reflectorSideCyl_solid,reflectorHollow_solid);
    reflectorSide_log = new G4LogicalVolume(reflectorFinal_solid,MaterialReflector,"ReflecSide_log",0,0,0);
  }

  
  // --- Teflon
  MaterialTeflon = G4Material::GetMaterial("G4_TEFLON");
  G4cout << "Inserting the Teflon " << G4endl 
         << "          section:  " << teflonSideSection 
         << "          height:   " << teflonSideHeight << G4endl;
  G4cout << "          material: " << MaterialTeflon->GetName() << G4endl;          


  if(type == "Box")
  {
    // --- Teflon side
    nameTemp = name + "TeflonS_sol";
    teflonSideBox_solid = new G4Box(nameTemp,teflonSideSection/2.0,teflonSideSection/2.0,teflonSideHeight/2.0);

    // The teflon side is a box with with a hollow obtained by subtraction of volumes
    // Teflon side hollow definition
    G4double teflonSideHollowSection = teflonSideSection - 2.0*teflonSideThickness;
    G4double teflonSideHollowHeight  = teflonSideHeight;
    G4Box *teflonHollow_solid = new G4Box("TeflonHollow_sol",teflonSideHollowSection/2.0,teflonSideHollowSection/2.0,teflonSideHollowHeight/2.0);
    
    // Now the final solid as a subtraction
    G4SubtractionSolid* teflonFinal_solid = new G4SubtractionSolid("TeflonSide_sol",teflonSideBox_solid,teflonHollow_solid);
    teflonSide_log = new G4LogicalVolume(teflonFinal_solid,MaterialTeflon,"TeflonSide_log",0,0,0);
    
    // Teflon front
    nameTemp = name + "TeflonFront_sol";
    teflonFrontBox_solid = new G4Box(nameTemp,teflonFrontSection/2.0,teflonFrontSection/2.0,teflonFrontHeight/2.0);
    teflonFront_log = new G4LogicalVolume(teflonFrontBox_solid,MaterialTeflon,"TeflonFront_log",0,0,0);
  }
  else  // assumed cylinder
  {
    // --- Teflon side
    nameTemp = name + "TeflonS_sol";
    G4double outerRadius  = teflonSideSection/2.0;
    G4double innerRadius  = 0.0;
    G4double initialAngle = 0.0*deg;
    G4double finalAngle   = 360.0*deg;
    teflonSideCyl_solid = new G4Tubs(nameTemp,innerRadius,outerRadius,teflonSideHeight/2.0,initialAngle,finalAngle);

    // Teflon side hollow definition
    outerRadius  = teflonSideSection/2.0 - teflonSideThickness;
    G4double teflonSideHollowHeight = teflonSideHeight;
    G4Tubs* teflonHollow_solid = new G4Tubs("teflonhollow",innerRadius,outerRadius,teflonSideHollowHeight/2.0,initialAngle,finalAngle);

    // Now the final solid as a subtraction
    G4SubtractionSolid* teflonFinal_solid = new G4SubtractionSolid("TeflonSide_sol",teflonSideCyl_solid,teflonHollow_solid);
    teflonSide_log = new G4LogicalVolume(teflonFinal_solid,MaterialTeflon,"TeflonSide_log",0,0,0);

    // Teflon front
    nameTemp = name + "TeflonFront_sol";
    outerRadius  = teflonSideSection/2.0;
    teflonFrontCyl_solid = new G4Tubs(nameTemp,innerRadius,outerRadius,teflonFrontHeight/2.0,initialAngle,finalAngle);
    teflonFront_log = new G4LogicalVolume(teflonFrontCyl_solid,MaterialTeflon,"TeflonFront_log",0,0,0);
  }
  

 // Placement of all parts in the Mother Volume ----------------------

 // --- Scintillator placement
 nameTemp = scintiName;
 new G4PVPlacement(0,scintillatorPosition,scintillator_log,scintiName,detMother_log,false,0);

 // --- Reflector placement
 nameTemp = name + "Reflector";
 new G4PVPlacement(0,reflectorSidePosition,reflectorSide_log,nameTemp,detMother_log,false,0);
 
 // --- Teflon placement
 // ------ side
 nameTemp = name + "TeflonSide";
 new G4PVPlacement(0,teflonSidePosition,teflonSide_log,nameTemp,detMother_log,false,0);
// ------ front
 nameTemp = name + "TeflonFront";
 new G4PVPlacement(0,teflonFrontPosition,teflonFront_log,nameTemp,detMother_log,false,0); 
 
 
// ----------------------------------------------- Color attributions
  scintillator_log->SetVisAttributes(blue_scintillator);
  //coverFront_log->SetVisAttributes(red_cover);
  //coverSide_log->SetVisAttributes(red_cover);
  teflonFront_log->SetVisAttributes(white_teflon);
  teflonSide_log->SetVisAttributes(white_teflon);
  reflectorSide_log->SetVisAttributes(yellow_reflector);
  
  G4cout << " Detector "<< name << " built" << G4endl;  
}

void DetectorElement::printInfo(std::ofstream *outputResults)
{
  *outputResults << G4endl;
  *outputResults << "--- Detector properties ---" << G4endl;
  *outputResults << "--- Detector Type:      " << type << G4endl;  
  *outputResults << "--- Scintillator ---" << G4endl;  
  *outputResults << " Scintillator material: " << MaterialScintillator->GetName() << G4endl;
 if(type == "Box")
 {
  *outputResults << " Section side:          " << scintillatorSection/mm << " mm" << G4endl;
 }
 else
 {
  *outputResults << " Section diameter:      " << scintillatorSection/mm << " mm" << G4endl;
 }
  *outputResults << " Height:                " << scintillatorHeight/mm << " mm" << G4endl;

  *outputResults << "--- Reflector ---" << G4endl;
  *outputResults << " Reflector material:    " << MaterialReflector->GetName() << G4endl;
  *outputResults << " Thickness:             " << reflectorSideThickness/mm << " mm" << G4endl;

  *outputResults << "--- Teflon ---" << G4endl;
  *outputResults << " Teflon material:       " << MaterialTeflon->GetName() << G4endl;
  *outputResults << " Side Thickness:        " << teflonSideThickness/mm << " mm" << G4endl;
  *outputResults << " Front Thickness:       " << teflonFrontHeight/mm << " mm" << G4endl;
  
/* 
  *outputResults << "--- Cover ---" << G4endl;
  *outputResults << " Cover material:        " << MaterialCover->GetName() << G4endl;
  *outputResults << " Side  thickness:       " << coverSideThickness/mm << " mm" << G4endl;
  *outputResults << " Front thickness:       " << coverFrontHeight/mm << " mm" << G4endl;
*/

  *outputResults << " --- Complete detector ---" << G4endl;
 if(type == "Box")
 {
  *outputResults << " Section side:          " << detMotherSection/mm << " mm" << G4endl;
 }
 else
 {
  *outputResults << " Section diameter:      " << detMotherSection/mm << " mm" << G4endl;    
 }
  *outputResults << " Height:                " << detMotherHeight/mm << " mm" << G4endl;
  *outputResults << G4endl;
  return;
}


  


