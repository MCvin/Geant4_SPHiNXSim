//
//  PhysicsList.hh 2013-09-10  Maxime Chauvin
//

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
   ~PhysicsList();

    virtual void SetCuts();
};

#endif
