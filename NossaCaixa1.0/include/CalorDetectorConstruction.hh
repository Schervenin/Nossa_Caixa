//*******************************************************
// CalorDetectorConstruction.hh


//*******************************************************

#ifndef CalorDetectorConstruction_h
#define CalorDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "DetectorElement.hh"

#include <map>

class G4Material;
class G4Box;
class G4Tubs;
class G4Sphere;
class G4LogicalVolume;
class G4VisAttributes;
class CalorDetectorMessenger;
class ofstream;

class CalorDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    CalorDetectorConstruction();
   ~CalorDetectorConstruction() override;

   G4VPhysicalVolume* Construct() override;
   
   CalorDetectorMessenger *detectorMessenger = nullptr;
      
  protected:
   void DefineMaterials();
   void DefineColors();
  
  protected:
   // World
   G4Box*             world_solid = nullptr;
   G4LogicalVolume*   world_logic = nullptr;
   G4VPhysicalVolume* world_phys  = nullptr;
   
   // detectors of the type DetectorElement
   G4double           detectorSection,detectorHeight;
   G4String           detectorMaterial, detectorShape;
   
   
   G4int              numDetectors, numCaps;
   G4String           detectorName;
   G4String           capName;
   G4ThreeVector      *detectorPosition = nullptr, *detectorAngle = nullptr;
   G4double           radialDislocation;  // to expand or contract the detectors radial distance
   G4String           scintiName;
   G4double           scintiSection,scintiHeight;
   G4ThreeVector      scintiPosition;  // to be used by PrimaryGenerator to produce the intrinsic radioactivity
   
//   G4double           *Xcm,*Ycm,*Zcm;
//   G4double           *Theta,*Phi,*Psi;
   DetectorElement    *detector = nullptr;
   
   std::map<G4int,G4String> DetectorIndex;  // map of (detector name,detector index)   
    std::map<G4int,G4String> CapIndex;  // map of (cap name,cap index)   
      
   G4String fileDirectory, detDescriptorName, fileName, capDescriptorName, fileNameCap;
   
  // The target holder 
   G4Material*        targetHolderMaterial  = nullptr;
   //G4Box*             boxTargetHolder_solid = nullptr;
   //G4double           widthBoxTH, lengthBoxTH, thicknessTH;  // X, Y, Z
   G4double           thicknessTH; 
   G4double           outerRadiusTH,innerRadiusTH;
   G4ThreeVector      targetHolderPosition;
   G4Tubs*            cylTargetHolder_solid = nullptr;
   G4LogicalVolume*   targetHolder_logic    = nullptr;
   G4VPhysicalVolume* targetHolder_phys     = nullptr;
   
  // Spherical part of the Chamber
   G4Material*         sphChamberMaterial = nullptr;
   G4double            sphChInnerRadius,sphChOuterRadius,sphChThickness,sphLength;
   G4double            AngleSphere;
   G4Tubs*           sphChamber_solid = nullptr;
   G4LogicalVolume*    sphChamber_logic = nullptr;
   G4ThreeVector       sphChamberPosition;      
   G4VPhysicalVolume*  sphChamber_phys;
   
  // prototype box for Cap (=> 0 <=)
   G4Material* CapMaterial = nullptr;
   G4SubtractionSolid*      Cap_solid = nullptr;
   G4double    widthCap, lengthCap, thicknessCap, wallxCap, wallyCap, frontCap;  // X, Y, Z
   G4ThreeVector    *capPosition = nullptr, *capAngle = nullptr;
   G4LogicalVolume*  Cap_logic    = nullptr;
   G4VPhysicalVolume* Cap_phys     = nullptr;
   std::vector<G4VPhysicalVolume*> CapsVolumes;

   
  public:
   DetectorElement* GetDetectorElement()     {return detector;}
   G4String GetDetectorShape()               {return detectorShape;}
   G4double GetDetectorSection()             {return detectorSection;}
   G4double GetDetectorHeight()              {return detectorHeight;}
   G4String GetScintillatorName()            {return scintiName;}
   G4double GetScintillatorSection()         {return scintiSection;}
   G4double GetScintillatorHeight()          {return scintiHeight;}
   G4ThreeVector GetScintillatorPosition()   {return scintiPosition;}
   
   
   G4int            GetNumberOfDetectors()   {return numDetectors;}
   G4int            GetNumberOfCaps()        {return numCaps;}
   G4String         GetNameOfDetectors()     {return detectorName;}
   
   G4ThreeVector*   GetDetectorPosition()    {return detectorPosition;}
   G4ThreeVector*   GetDetectorAngle()       {return detectorAngle;}
   
   G4String GetDetectorsFileDirectory()      {return fileDirectory;}
   G4String GetDetectorsDescriptorName()     {return detDescriptorName;}
   G4String GetDetectorsFileName()           {return fileName;}
   
      
   void     SetDetectorsFileDirectory(G4String dir);  
   void     SetDetectorsDescriptorName(G4String name);
   void     SetDetectorsMaterial(G4String name);
   void     SetDetectorsShape(G4String name);
   void     SetRadialDislocation(G4double dR);  // expansion (+) or contraction (-) for detector position

   void     SetChamberMaterial(G4String name);   
   void     SetChamberThickness(G4double thick);
   void     SetChamberOuterRadius(G4double outRadius);
   void     SetTargetHolderMaterial(G4String name);
   

   std::map<G4int,G4String>*  GetDetectorsMapPointer() {return &DetectorIndex;}
   
   void  printDetectorInfo(std::ofstream *);

  private:
   // Colors
   G4VisAttributes* grey_wire  = nullptr;   // for vacuum
   G4VisAttributes* red_copper = nullptr;  // for the target support
   G4VisAttributes* grey_solid = nullptr;  // for the spherical part of the chamber
   
   void insertDetectors();
   void insertCaps();
};
#endif

