//
//  DetectorMessenger.cc 2013-12-18  Maxime Chauvin
//

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* det): detector(det)
{
  YRotationCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetInstRotY",this);
  YRotationCmd->SetGuidance("Set Y rotation of the instrument");
  YRotationCmd->SetParameterName("Rotation Y",false);
  YRotationCmd->SetUnitCategory("Angle");
  YRotationCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetFastOuterCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetFastOuter",this);
  SetFastOuterCmd->SetGuidance("Set plastic diameter");
  SetFastOuterCmd->SetParameterName("Plastic diameter",false);
  SetFastOuterCmd->SetUnitCategory("Length");
  SetFastOuterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetFastHeightCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetFastHeight",this);
  SetFastHeightCmd->SetGuidance("Set plastic height");
  SetFastHeightCmd->SetParameterName("Plastic height",false);
  SetFastHeightCmd->SetUnitCategory("Length");
  SetFastHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetHighZThickCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetHighZThick",this);
  SetHighZThickCmd->SetGuidance("Set high Z thickness");
  SetHighZThickCmd->SetParameterName("HighZ thickness",false);
  SetHighZThickCmd->SetUnitCategory("Length");
  SetHighZThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetHighZHeightCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetHighZHeight",this);
  SetHighZHeightCmd->SetGuidance("Set high Z height");
  SetHighZHeightCmd->SetParameterName("HighZ height",false);
  SetHighZHeightCmd->SetUnitCategory("Length");
  SetHighZHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetShieldHeightCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetShieldHeight",this);
  SetShieldHeightCmd->SetGuidance("Set shield height");
  SetShieldHeightCmd->SetParameterName("Shield height",false);
  SetShieldHeightCmd->SetUnitCategory("Length");
  SetShieldHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetShieldThick_PbCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetShieldThick_Pb",this);
  SetShieldThick_PbCmd->SetGuidance("Set shield Pb thickness");
  SetShieldThick_PbCmd->SetParameterName("SetShieldThick_Pb thickness",false);
  SetShieldThick_PbCmd->SetUnitCategory("Length");
  SetShieldThick_PbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetShieldThick_SnCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetShieldThick_Sn",this);
  SetShieldThick_SnCmd->SetGuidance("Set shield Sn thickness");
  SetShieldThick_SnCmd->SetParameterName("SetShieldThick_Sn thickness",false);
  SetShieldThick_SnCmd->SetUnitCategory("Length");
  SetShieldThick_SnCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetShieldThick_CuCmd = new G4UIcmdWithADoubleAndUnit("/detector/SetShieldThick_Cu",this);
  SetShieldThick_CuCmd->SetGuidance("Set shield Cu thickness");
  SetShieldThick_CuCmd->SetParameterName("SetShieldThick_Cu thickness",false);
  SetShieldThick_CuCmd->SetUnitCategory("Length");
  SetShieldThick_CuCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete YRotationCmd;
  delete SetFastOuterCmd;
  delete SetFastHeightCmd;
  delete SetHighZThickCmd;
  delete SetHighZHeightCmd;
  delete SetShieldHeightCmd;
  delete SetShieldThick_PbCmd;
  delete SetShieldThick_SnCmd;
  delete SetShieldThick_CuCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == YRotationCmd ) detector->SetInstRotY(YRotationCmd->GetNewDoubleValue(newValue));
  if( command == SetFastOuterCmd ) detector->SetFastOuter(SetFastOuterCmd->GetNewDoubleValue(newValue));
  if( command == SetFastHeightCmd ) detector->SetFastHeight(SetFastHeightCmd->GetNewDoubleValue(newValue));
  if( command == SetHighZThickCmd ) detector->SetHighZThick(SetHighZThickCmd->GetNewDoubleValue(newValue));
  if( command == SetHighZHeightCmd ) detector->SetHighZHeight(SetHighZHeightCmd->GetNewDoubleValue(newValue));
  if( command == SetShieldHeightCmd ) detector->SetShieldHeight(SetShieldHeightCmd->GetNewDoubleValue(newValue));
  if( command == SetShieldThick_PbCmd ) detector->SetShieldThick_Pb(SetShieldThick_PbCmd->GetNewDoubleValue(newValue));
  if( command == SetShieldThick_SnCmd ) detector->SetShieldThick_Sn(SetShieldThick_SnCmd->GetNewDoubleValue(newValue));
  if( command == SetShieldThick_CuCmd ) detector->SetShieldThick_Cu(SetShieldThick_CuCmd->GetNewDoubleValue(newValue));
}
