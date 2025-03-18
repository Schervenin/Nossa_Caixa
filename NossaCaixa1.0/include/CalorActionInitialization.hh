// ActionInitialization.hh
//

#ifndef CalorActionInitialization_h
#define CalorActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
//#include "CalorDetectorConstruction.hh"
//#include "CalorPrimaryGeneratorAction.hh"

class CalorDetectorConstruction;
class CalorPrimaryGeneratorAction;

//________________________________________________________________________
class CalorActionInitialization : public G4VUserActionInitialization
{
  public:
    CalorActionInitialization(CalorDetectorConstruction* );
    virtual ~CalorActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    CalorDetectorConstruction*    DetCon;
//    CalorPrimaryGeneratorAction*  PrimGen;
};

#endif

    
