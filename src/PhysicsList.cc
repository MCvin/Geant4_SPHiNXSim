//
//  PhysicsList.cc 2013-09-10  Maxime Chauvin
//

#include "PhysicsList.hh"

#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  // EM Livermore Polarized physics
  RegisterPhysics(new G4EmLivermorePolarizedPhysics());

  defaultCutValue = 0.2*mm; //10% of the slow scintillator thickness (2mm)
  SetVerboseLevel(1);
}

PhysicsList::~PhysicsList(){}

void PhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
  // fix lower limit for cut
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV, 1*GeV);

  // print the table of elements with the cuts
  DumpCutValuesTable();
}  
