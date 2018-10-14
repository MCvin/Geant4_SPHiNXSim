//
//  DetectorMessenger.hh 2013-12-18  Maxime Chauvin
//

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

class DetectorConstruction;

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    DetectorConstruction* detector;

    G4UIcmdWithADoubleAndUnit* YRotationCmd;
    G4UIcmdWithADoubleAndUnit* SetFastOuterCmd;
    G4UIcmdWithADoubleAndUnit* SetFastHeightCmd;
    G4UIcmdWithADoubleAndUnit* SetHighZThickCmd;
    G4UIcmdWithADoubleAndUnit* SetHighZHeightCmd;
    G4UIcmdWithADoubleAndUnit* SetShieldHeightCmd;
    G4UIcmdWithADoubleAndUnit* SetShieldThick_PbCmd;
    G4UIcmdWithADoubleAndUnit* SetShieldThick_SnCmd;
    G4UIcmdWithADoubleAndUnit* SetShieldThick_CuCmd;
};

#endif

