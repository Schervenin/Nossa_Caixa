//  Calor Detector Construction

// version 3.0 - 26/06/2019
//   Inclusion of the option of cylindrical geometry for the detectors
//   Inclusion of a target holder as described by Calvo in the meeting of April 2019
//   Inclusion of the spherical part of the scattering chamber  

#include "CalorDetectorConstruction.hh"
#include "CalorDetectorMessenger.hh"
#include "DetectorElement.hh" 

#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Element.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"

//#include "G4BREPSolid.hh"
//#include "G4Polygon.hh"
//#include "G4Surface.hh"

#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"   

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


#include <fstream>     // to save the report file
#include <cmath>
#include <map>
#include <vector>

//___________________________________________________________________________
// Functions
// Conversion of number to string
template <typename T>
G4String NumberToString ( T Number )
{
	std::stringstream ss;
	ss << Number;
	return ss.str();
}

//___________________________________________________________________________
CalorDetectorConstruction::CalorDetectorConstruction()
{
  G4cout << "Building the DetectorConstruction..." << G4endl;
  DefineMaterials(); // before all
  DefineColors();

  // Target holder initial parameters
  targetHolderMaterial = G4Material::GetMaterial("G4_Cu");  // copper can be changed with macro command
  //widthBoxTH           = 20.0*mm;
  //lengthBoxTH          = 40.0*mm;
  thicknessTH          = 0.5*mm;
  innerRadiusTH        = 3.0*mm;
  outerRadiusTH        = 8.0*mm;
  targetHolderPosition.set(0.0,0.0,0.0);  // target at origin
  
  // Cap initial parameters
  CapMaterial = G4Material::GetMaterial("G4_Al");   // aluminum
  widthCap          = 43.9*mm;
  lengthCap         = 49.9*mm;
  thicknessCap      = 64.0*mm;
  wallxCap          = 2.0*mm;
  wallyCap          = 5.0*mm;
  frontCap          = 4.0*mm;
  //capPosition.set(0.0,0.0,82.8);  // somewhere
  
  // Spherical top of the scattering chamber
  sphChamberMaterial = G4Material::GetMaterial("G4_Al");  // aluminum
  sphChOuterRadius   = 200.0*mm;
  sphChThickness     = 5.0*mm;
  sphLength          = 10.*cm;
  sphChInnerRadius   = sphChOuterRadius - sphChThickness;
  AngleSphere        = 90.*degree;  // angle theta covered by the spherical top
  sphChamberPosition.set(0.0,0.0,0.0);  // center of the sphere at the target position
  
  detectorName = "D_";
  capName = "Cap_";
  // Default values
  detectorMaterial       = "LYSO";      // other option: LaBr3, GAGG (can be changed with macro command)
  detectorShape          = "Box";       // other option: Cylinder (can be changed with macro command)
  
  detector = new DetectorElement(detectorMaterial,detectorShape,detectorName);
  detectorHeight    = detector->getDetectorHeight();
//  radialDislocation = detectorHeight/2.0;
  radialDislocation = 0.0;  
  scintiName = detector->getScintillatorName();
  
  // Choose the file that contains rows with Phi, Theta, Psi of the detectors
  // The first line of the file must contain the number of detectors
  fileDirectory      = "../DetectorsPositionFiles/";         
  detDescriptorName  = "lplaceNC.dat";    // NC detectors
//  capDescriptorName  = "lplaceNCcaps.dat";    // NC caps
  capDescriptorName  = "lplaceCaps.dat";    // NC caps provisório
//  detDescriptorName      = "lplace-radial.dat";      // radial configuration
  
  detectorMessenger = new CalorDetectorMessenger(this);
}

//___________________________________________________________________________
CalorDetectorConstruction::~CalorDetectorConstruction() 
{
  delete detectorMessenger;
}

//___________________________________________________________________________
G4VPhysicalVolume* CalorDetectorConstruction::Construct()
{
  G4cout << G4endl << "DetectorConstruction::Construct()" << G4endl;

  //----- Build world
  G4double worldXDimension = 3.0*m;
  G4double worldYDimension = 3.0*m;
  G4double worldZDimension = 3.0*m;

  G4Material *material;
  material = G4Material::GetMaterial("G4_Galactic");  // vacuum
  if(world_solid)
  {
    delete world_solid;
    world_solid = 0;
  }
  if(world_logic)
  {
    delete world_logic;
    world_logic = 0;
  }

  world_solid = new G4Box( "WorldSolid",worldXDimension/2.0,worldYDimension/2.0,worldZDimension/2.0 );
  world_logic = new G4LogicalVolume( world_solid, material, "WorldLogical", 0, 0, 0 );
  if(world_phys) 
  {
    delete detector;
    detector = 0;
    delete world_phys;
    detector = new DetectorElement(detectorMaterial,detectorShape,detectorName);
  }
  world_phys  = new G4PVPlacement( 0,G4ThreeVector(0,0,0),"World",world_logic,0,false,0);
  
  // --- insert prototype Cap (=> 0 <=)
  if(Cap_solid)
  {
	  delete Cap_solid;
	  Cap_solid=0;
  }
  
  
 /* 
 //**************** <=0=> 
    G4Box *CapBulk_solid = new G4Box( "CapBulk",widthCap/2.0,lengthCap/2.0,thicknessCap/2.0 );
 

  G4Tubs* cyl1= new G4Tubs("Cylinder1",0*mm,12*mm,10*mm,0,360.*deg); 
 
  G4ThreeVector relativeCutPosition(0.0 * mm, 0.0 * mm, -35*mm); // Translation
  G4RotationMatrix* rot1 = new G4RotationMatrix();
  rot1->rotateX(45.0*deg);
  rot1->rotateY(45.0*deg);
  G4SubtractionSolid * cut1 = new G4SubtractionSolid("cut1",CapBulk_solid,cyl1,rot1,relativeCutPosition);

  
  G4Box *CapHollow_solid = new G4Box( "CapBulk", widthCap/2.0 - wallxCap, lengthCap/2.0 - wallyCap, thicknessCap/2.0 - frontCap);
  G4ThreeVector relativeHollowPosition(0.0 * mm, 0.0 * mm, frontCap/2.); // Translation
  //Cap_solid = new G4SubtractionSolid("Cap_sol",CapBulk_solid,CapHollow_solid,0,relativeHollowPosition);

  Cap_solid = new G4SubtractionSolid("Cap_sol",cut1,CapHollow_solid,0,relativeHollowPosition);
 //**************** <=0=>  
  */
  
  auto mesh = CADMesh::TessellatedMesh::FromSTL("SoDetBox.stl");
  Cap_solid = mesh->GetSolid();
    
  if(Cap_logic)
  {
    delete Cap_logic;
    Cap_logic = 0;
  }
  Cap_logic = new G4LogicalVolume(Cap_solid,CapMaterial,"Cap_log",0,0,0);

  //Cap_phys = new G4PVPlacement(0,G4ThreeVector(0,0,50.0*mm),Cap_logic,"Cap",world_logic,false,0);

  /* / --- Insert the Target Holder
  if(boxTargetHolder_solid)
  {
    delete boxTargetHolder_solid;
    boxTargetHolder_solid = 0;
  }
  if(cylTargetHolder_solid)
  {
    delete cylTargetHolder_solid;
    cylTargetHolder_solid = 0;
  }*/
  
  //boxTargetHolder_solid = new G4Box( "BoxTHSolid",widthBoxTH/2.0,lengthBoxTH/2.0,thicknessTH/2.0 );
  G4double initialAngle = 0.0*deg;
  G4double finalAngle   = 360.0*deg;
  G4Tubs* targetHolder_solid = new G4Tubs("CylTHSolid",innerRadiusTH,outerRadiusTH,thicknessTH/2.0,initialAngle,finalAngle);
  // The target holder is the Union of these 2 solids, the Box displaced by some distace
  //G4double boxDisplacement = lengthBoxTH/2.0+outerRadiusTH - 5.0*mm;
  //G4ThreeVector zTrans(0, boxDisplacement, 0);
  //G4UnionSolid* targetHolder_solid = new G4UnionSolid("Cyl+BoxMoved", cylTargetHolder_solid, boxTargetHolder_solid, 0, zTrans);
  
  if(targetHolder_logic)
  {
    delete targetHolder_logic;
    targetHolder_logic = 0;
  }
  targetHolder_logic = new G4LogicalVolume(targetHolder_solid,targetHolderMaterial,"targetHolder_log",0,0,0);
  targetHolder_phys = new G4PVPlacement(0,targetHolderPosition,targetHolder_logic,"targetHolder",world_logic,false,0);
  
  // --- Insert the spherical top of the scattering chamber
  G4double PhiInitial   = 0.0*degree;
  G4double DeltaPhi     = 360.0*degree;
  G4double ThetaInitial = 0.0*degree;
  G4double DeltaTheta   = AngleSphere;

  if(sphChamber_solid) 
  {  
    delete sphChamber_solid;
    sphChamber_solid = 0;
  }
//  sphChamber_solid = new G4Sphere("sphChamber_solid",sphChInnerRadius,sphChOuterRadius,PhiInitial,DeltaPhi,ThetaInitial,DeltaTheta);
 // sphChamber_solid = new G4Tubs("sphChamber_solid",sphChInnerRadius,sphChOuterRadius,sphLength,0.,360.);
  
  auto meshSup = CADMesh::TessellatedMesh::FromSTL("SoSuportes.stl");
  sphChamber_solid = meshSup->GetSolid();
  

  if(sphChamber_logic) 
  {  
    delete sphChamber_logic;
    sphChamber_logic = 0;
  }
  sphChamber_logic = new G4LogicalVolume(sphChamber_solid,sphChamberMaterial,"sphChamber_log",0,0,0);

  G4RotationMatrix* xrot = new G4RotationMatrix();
  xrot->rotateX(90.0*deg);
    
  sphChamber_phys  = new G4PVPlacement(xrot,sphChamberPosition,sphChamber_logic,"sphChamber",world_logic,false,0);

  
// --- World invisible
  world_logic->SetVisAttributes(G4VisAttributes::GetInvisible);  // invisible
//  world_logic->SetVisAttributes(grey_wire);
  if(targetHolderMaterial->GetName() == "G4_Cu")
    targetHolder_logic->SetVisAttributes(red_copper);
  else if(targetHolderMaterial->GetName() == "G4_Galactic")
      targetHolder_logic->SetVisAttributes(grey_wire);
  else
    targetHolder_logic->SetVisAttributes(grey_solid);
  
  if(sphChamberMaterial->GetName() == "G4_Galactic")
      sphChamber_logic->SetVisAttributes(grey_wire);
  else
      sphChamber_logic->SetVisAttributes(grey_solid);

  // Choose the detector
  if(detector)  // to apply changes
  {
    delete detector;  
    detector = 0;
  }
  detector = new DetectorElement(detectorMaterial,detectorShape,detectorName);
  scintiSection = detector->getScintillatorSection();
  scintiHeight  = detector->getScintillatorHeight();

  insertDetectors();
  insertCaps();
  
  G4cout << G4endl;
  
  return world_phys;
}


//___________________________________________________________________________
void CalorDetectorConstruction::insertDetectors()
{
  G4cout << " Inserting Detectors... " << G4endl;
  
  if(detectorPosition)  // to apply changes
  {
    delete[] detectorPosition;  
    detectorPosition = 0;
  }
  if(detectorAngle)  // to apply changes
  {
    delete[] detectorAngle;  
    detectorAngle = 0;
  }

  // Choose the file that contains rows with Theta, Phi, Psi of the detectors
  // The first line of the file must contain the number of detectors
  
  fileName = fileDirectory + detDescriptorName;
  G4cout << " File for the calorimeter's detectors positions:" << G4endl;
  G4cout << fileName << G4endl;
  
  std::ifstream detDescriptorFile;
  detDescriptorFile.open(fileName,std::ios::in);
  if(!detDescriptorFile)
  {
      G4cout << "Can not open file " << fileName << G4endl;
      G4cout << "Only one detector will be inserted on the z axis" << G4endl;
      G4cout << "at +30 cm from the target/source position" << G4endl;
      
      detectorPosition    = new G4ThreeVector[1];
      detectorAngle       = new G4ThreeVector[1];
      detectorPosition[0] = G4ThreeVector(0.,0.,300.0*mm);
      detectorAngle[0]    = G4ThreeVector(0.,0.,0.);
      
      G4String positionName = NumberToString(0)+"_"+NumberToString(0);
      G4String thisCopyName = detectorName+positionName;
 
    new G4PVPlacement(0,detectorPosition[0], 
				      detector->getDetectorPointer(),thisCopyName,
				      world_logic,false,0);
     DetectorIndex[0] = thisCopyName;            // adds the index in the map
     G4cout << "Detector Name: " << thisCopyName << " Index: " << 0 << G4endl;     
  }
  else
  {
      detDescriptorFile >> numDetectors;  // Reads the first entry of the file: numbers of detectors
      G4cout << " Numbers of detectors for this configuration: " << numDetectors << G4endl;
      if(numDetectors <= 0)
      {
          G4cout << " Problem with the detector file " << fileName << G4endl;
          G4cout << " No detectors to be inserted." << G4endl;
          return;
      }
      detectorPosition  = new G4ThreeVector[numDetectors];
      detectorAngle     = new G4ThreeVector[numDetectors];
          
      G4cout << "__________________________________________________________________________" << G4endl;
      G4cout << "Det   X        Y        Z        Phi      Theta    Psi        Det. Name" << G4endl;
      // reads the rest of the file: numDetectors rows with X, Y, Z, Theta, Phi, Psi
      for (G4int copyNr = 0; copyNr < numDetectors; copyNr++)  
      {
        G4double posX,posY,posZ,Phi,Theta,Psi;
        detDescriptorFile >>  posX >> posY >> posZ >> Phi >> Theta >> Psi;
        
        // For dislocation in the radial diraction (expansion of the ball)
        // To avoid overlaps of detector elements
        G4ThreeVector originalPosition = G4ThreeVector(posX,posY,posZ); // as read from the position file
        G4double      norm = sqrt(posX*posX+posY*posY+posZ*posZ); 
        G4ThreeVector radialVersor = originalPosition/norm;       // the radial versor
 
        G4ThreeVector newPosition = originalPosition+radialDislocation*radialVersor;
        
        // To test volume overlaping: /geometry/test/run
        
        // new position coordinates
        posX = newPosition.x();  
        posY = newPosition.y();
        posZ = newPosition.z();

        G4cout.flags(std::ios::left);
        G4cout.setf(std::ios::floatfield,std::ios::fixed);
        G4cout << std::setfill(' ') 
                << std::setw(5) << copyNr                << " "
                << std::setprecision(2)
                << std::setw(8) << posX/mm        << " "
                << std::setw(8) << posY/mm        << " "
                << std::setw(8) << posZ/mm        << " "
                << std::setw(8) << Phi/degree    << " "          
                << std::setw(8) << Theta/degree  << " "
                << std::setw(8) << Psi/degree    << "   ";

        detectorPosition[copyNr] = G4ThreeVector(posX,posY,posZ);
        detectorAngle[copyNr]    = G4ThreeVector(Phi,Theta,Psi);

     //rotationMat see example http://geant4.web.cern.ch/geant4/UserDocumentation/Doxygen/examples_doc/html_transforms/html/DetectorConstruction_8cc_source.html

         // The line below works for the new Zero's file format
         G4RotationMatrix rotm1Inv = G4RotationMatrix(Phi+pi/2.0,Theta,Psi);       
         G4RotationMatrix rotm1 = rotm1Inv.inverse();

         G4Transform3D transform1 = G4Transform3D(rotm1,detectorPosition[copyNr]); // rotation + translation
          
         G4String positionName = NumberToString(round(Phi/degree))+"_"+NumberToString(round(Theta/degree));
         G4String thisCopyName = detectorName+NumberToString(copyNr)+"_"+positionName;
 
         new G4PVPlacement(transform1, 
                          detector->getDetectorPointer(),thisCopyName,
                          world_logic,false,copyNr);
         DetectorIndex[copyNr] = thisCopyName;  // adds the index in the map
         G4cout << thisCopyName << G4endl;     
      }
      detDescriptorFile.close();
      
      G4cout << "_________________________________________________________________________" << G4endl;
      G4cout << "    ... Detectors inserted. " << G4endl;
  }
  
  
}

//___________________________________________________________________________
void CalorDetectorConstruction::insertCaps()
{
  G4cout << " Inserting Caps... " << G4endl;
  
  if(capPosition)  // to apply changes
  {
    delete[] capPosition;  
    capPosition = 0;
  }
  if(capAngle)  // to apply changes
  {
    delete[] capAngle;  
    capAngle = 0;
  }

  // Choose the file that contains rows with Theta, Phi, Psi of the detectors
  // The first line of the file must contain the number of detectors
  
  fileNameCap = fileDirectory + capDescriptorName;
  G4cout << " File for the NC detector cap positions:" << G4endl;
  G4cout << fileNameCap << G4endl;
  
  std::ifstream capDescriptorFile;
  capDescriptorFile.open(fileNameCap,std::ios::in);
  if(!capDescriptorFile)
  {
      G4cout << "Can not open file " << fileNameCap << G4endl;
      G4cout << "Only one cap will be inserted on the z axis" << G4endl;
      G4cout << "at +30 cm from the target/source position" << G4endl;
      
    capPosition    = new G4ThreeVector[1];
    capAngle       = new G4ThreeVector[1];
      capPosition[0]=  G4ThreeVector(0.,0.,300.0*mm);
      capAngle[0]= G4ThreeVector(0.,0.,0.);
      
     Cap_phys = new G4PVPlacement(0,capPosition[0],Cap_logic,"Cap",world_logic,false,0);
	//Cap_phys = new G4PVPlacement(0,G4ThreeVector(0,0,50.0*mm),Cap_logic,"Cap",world_logic,false,0);     
	
     // G4String positionName = NumberToString(0)+"_"+NumberToString(0);
     // G4String thisCopyName = capName+positionName;
   
  }
  else
  {
      capDescriptorFile >> numCaps;  // Reads the first entry of the file: numbers of detectors
      G4cout << " Numbers of caps for this configuration: " << numCaps << G4endl;
      if(numCaps <= 0)
      {
          G4cout << " Problem with the cap file " << fileNameCap << G4endl;
          G4cout << " No caps to be inserted." << G4endl;
          return;
      }
      capPosition  = new G4ThreeVector[numCaps];
      capAngle     = new G4ThreeVector[numCaps];
          
      G4cout << "__________________________________________________________________________" << G4endl;
      G4cout << "Cap   X        Y        Z        Phi      Theta    Psi        Det. Name" << G4endl;
      // reads the rest of the file: numDetectors rows with X, Y, Z, Theta, Phi, Psi
      for (G4int copyNr = 0; copyNr < numCaps; copyNr++)  
      {
        G4double posX,posY,posZ,Phi,Theta,Psi;
        capDescriptorFile >>  posX >> posY >> posZ >> Phi >> Theta >> Psi;
        
        // For dislocation in the radial diraction (expansion of the ball)
       // To avoid overlaps of detector elements
        G4ThreeVector originalPosition = G4ThreeVector(posX,posY,posZ); // as read from the position file
        G4double      norm = sqrt(posX*posX+posY*posY+posZ*posZ); 
        G4ThreeVector radialVersor = originalPosition/norm;       // the radial versor
 
        G4ThreeVector newPosition = originalPosition+radialDislocation*radialVersor;
        
        // To test volume overlaping: /geometry/test/run
        
        // new position coordinates
        posX = newPosition.x();  
        posY = newPosition.y();
        posZ = newPosition.z();

        G4cout.flags(std::ios::left);
        G4cout.setf(std::ios::floatfield,std::ios::fixed);
        G4cout << std::setfill(' ') 
                << std::setw(5) << copyNr                << " "
                << std::setprecision(2)
                << std::setw(8) << posX/mm        << " "
                << std::setw(8) << posY/mm        << " "
                << std::setw(8) << posZ/mm        << " "
                << std::setw(8) << Phi/degree    << " "          
                << std::setw(8) << Theta/degree  << " "
                << std::setw(8) << Psi/degree    << "   ";

        capPosition[copyNr] = G4ThreeVector(posX,posY,posZ);
        capAngle[copyNr]    = G4ThreeVector(Phi,Theta,Psi);

     //rotationMat see example http://geant4.web.cern.ch/geant4/UserDocumentation/Doxygen/examples_doc/html_transforms/html/DetectorConstruction_8cc_source.html

         // The line below works for the new Zero's file format
         G4RotationMatrix rotm1Inv = G4RotationMatrix(Phi+pi/2.0,Theta,Psi);       
         G4RotationMatrix rotm1 = rotm1Inv.inverse();
         G4Transform3D transform1 = G4Transform3D(rotm1,capPosition[copyNr]); // rotation + translation
          
         G4String positionName = NumberToString(round(Phi/degree))+"_"+NumberToString(round(Theta/degree));
         G4String thisCopyName = capName+NumberToString(copyNr)+"_"+positionName;
         
		  new G4PVPlacement(transform1,Cap_logic,thisCopyName,world_logic,false,copyNr);
		  CapsVolumes.push_back(Cap_phys);
		  
		 // --> Cap_phys = new G4PVPlacement(0,G4ThreeVector(0,0,50.0*mm),Cap_logic,"Cap",world_logic,false,0);
         CapIndex[copyNr] = thisCopyName;  // adds the index in the map
         G4cout << thisCopyName << G4endl;     
      }
      capDescriptorFile.close();
      
      G4cout << "_________________________________________________________________________" << G4endl;
      G4cout << "    ... Caps inserted. " << G4endl;
  }
  
  
} 

//___________________________________________________________________________
void CalorDetectorConstruction::SetDetectorsFileDirectory(G4String dir)  
{
    G4cout << "The file directory for detectors description file was changed." << G4endl;
    G4cout << "Old file directory: " << fileDirectory << G4endl;
    fileDirectory = dir;
    G4cout << "New file directory: " << fileDirectory << G4endl;    
    fileName = fileDirectory + detDescriptorName;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetDetectorsDescriptorName(G4String name)  
{
    G4cout << "The file name for detectors was changed." << G4endl;
    G4cout << "Old file name: " << detDescriptorName << G4endl;
    detDescriptorName = name;
    G4cout << "New file name: " << detDescriptorName << G4endl;
    fileName = fileDirectory + detDescriptorName;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetDetectorsMaterial(G4String mat)  
{
    G4cout << "The detectors Material was changed." << G4endl;
    G4cout << "Old Material: " << detectorMaterial << G4endl;
    detectorMaterial = mat;
    G4cout << "New Material: " << detectorMaterial << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->PhysicsHasBeenModified();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetDetectorsShape(G4String shape)  
{
    G4cout << "The detectors Shape was changed." << G4endl;
    G4cout << "Old Shape: " << detectorShape << G4endl;
    detectorShape = shape;
    G4cout << "New Shape: " << detectorShape << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetRadialDislocation(G4double dR)  
{
    G4cout << "The detectors Radial Dislocation was changed." << G4endl;
    G4cout << "Old Radial Dislocation: " << radialDislocation/mm << " mm" << G4endl;
    radialDislocation = dR;
    G4cout << "New Radial Dislocation: " << radialDislocation/mm << " mm" << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetChamberMaterial(G4String mat)  
{
    G4cout << "The Chamber Material was changed." << G4endl;
    G4cout << "Old Material: " << sphChamberMaterial->GetName() << G4endl;
    sphChamberMaterial = G4Material::GetMaterial(mat);
    G4cout << "New Material: " << sphChamberMaterial->GetName() << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->PhysicsHasBeenModified();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetChamberThickness(G4double thick)  
{
    G4cout << "The chamber Thickness was changed." << G4endl;
    G4cout << "Old Thickness: " << sphChThickness/mm << " mm" << G4endl;
    sphChThickness = thick;
    G4cout << "New Thickness: " << sphChThickness/mm << " mm" << G4endl;
    sphChInnerRadius   = sphChOuterRadius - sphChThickness;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetChamberOuterRadius(G4double outRadius)  
{
    G4cout << "The chamber External Radius was changed." << G4endl;
    G4cout << "Old Radius: " << sphChOuterRadius/mm << " mm" << G4endl;
    sphChOuterRadius = outRadius;
    G4cout << "New Radius: " << sphChOuterRadius/mm << " mm" << G4endl;
    sphChInnerRadius = sphChOuterRadius - sphChThickness;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::SetTargetHolderMaterial(G4String mat)  
{
    G4cout << "The target holder Material was changed." << G4endl;
    G4cout << "Old Material: " << targetHolderMaterial->GetName() << G4endl;
    targetHolderMaterial = G4Material::GetMaterial(mat);
    G4cout << "New Material: " << targetHolderMaterial->GetName() << G4endl;
    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->PhysicsHasBeenModified();
    runManager->ReinitializeGeometry();
}

//___________________________________________________________________________
void CalorDetectorConstruction::DefineMaterials()
{
    G4cout << "DetectorConstruction: DefineMaterials." << G4endl;
    G4NistManager* NistMan = G4NistManager::Instance();
    G4bool isotopes = false;
    if(!G4Material::GetMaterial("G4_AIR",false))
      NistMan->FindOrBuildMaterial("G4_AIR", isotopes);
    if(!G4Material::GetMaterial("G4_Galactic",false))
      NistMan->FindOrBuildMaterial("G4_Galactic",isotopes);
    if(!G4Material::GetMaterial("G4_Cu",false))
      NistMan->FindOrBuildMaterial("G4_Cu",isotopes);
    if(!G4Material::GetMaterial("G4_Al",false))
      NistMan->FindOrBuildMaterial("G4_Al",isotopes);
}

//___________________________________________________________________________
void CalorDetectorConstruction::DefineColors()
{
    G4cout << "DetectorConstruction: DefineColors." << G4endl;

    grey_wire =  new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.2)); // 0,2 = transparência 
    grey_wire->SetVisibility(true);
    grey_wire->SetForceWireframe(true);    // gray wireframe
    grey_wire->SetLineWidth(1.0);
    
    red_copper   = new G4VisAttributes(G4Colour(1.0,0.4,0.4,1.0));
    red_copper->SetVisibility(true);
    red_copper->SetForceSolid(true);       // red solid
    
    grey_solid =  new G4VisAttributes(G4Colour(0.8,0.8,0.8,1.0)); // 0,5 = transparência 
    grey_solid->SetVisibility(true);
    grey_solid->SetForceSolid(true);      // gray solid
}

//___________________________________________________________________________
void CalorDetectorConstruction::printDetectorInfo(std::ofstream *output)
{
    *output << G4endl << " --- Chamber configuration" << G4endl;
    *output << " Material:     " << sphChamberMaterial->GetName() << G4endl;
    *output << " Outer Radius: " << sphChOuterRadius/mm << " mm" << G4endl;
    *output << " Thickness:    " << sphChThickness/mm << " mm" << G4endl;
    
    *output << G4endl << " --- Target holder" << G4endl;
    *output << " Material: " << targetHolderMaterial->GetName() << G4endl;
    
  
    *output << G4endl << " --- Calorimeter configuration" << G4endl;
    detector->printInfo(output);
    *output << " Directory of file for detectors positioning: " << G4endl << " " 
            << fileDirectory << G4endl;
    *output << " File with the detectors positioning:  " << detDescriptorName << G4endl;
    *output << " Calorimeter number of detectors:      " << numDetectors << G4endl;
    if(numDetectors == 1)
      *output << " Detector position: " << detectorPosition[0]/cm << " cm" << G4endl; 
    *output << " Radial extra dislocation:             " << radialDislocation << G4endl;
}

//___________________________________________________________________________


