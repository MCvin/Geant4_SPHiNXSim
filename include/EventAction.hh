//
//  EventAction.hh 2013-09-04  Maxime Chauvin
//

#ifndef EventAction_h
#define EventAction_h 1

#include "TrackerHit.hh"

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
    void CreateASCIIFile(const G4String name) {asciiFile = new std::ofstream(name);}
    void CreateROOTFile(const G4String name);

  private:
    struct Prim_event{
      int Id;
      float E;
      float Px, Py, Pz;
      float Dx, Dy, Dz;
    };
    struct Det_event{
	  int nhits; //Number of hits determines the size of arrays when writing the TTree
	  float Edep[1000];
	  float Px[1000];
	  float Py[1000];
	  float Pz[1000];
	  int UnitId[1000];
	  int nUnits;
	  float EdepTot[60];
      inline void SetAddress(TBranch* branch){
        TObjArray* leaf_array = branch->GetListOfLeaves();
        ((TLeaf*)leaf_array->At(0))->SetAddress(&nhits);
        ((TLeaf*)leaf_array->At(1))->SetAddress(Edep);
        ((TLeaf*)leaf_array->At(2))->SetAddress(Px);
        ((TLeaf*)leaf_array->At(3))->SetAddress(Py);
        ((TLeaf*)leaf_array->At(4))->SetAddress(Pz);
        ((TLeaf*)leaf_array->At(5))->SetAddress(UnitId);
        ((TLeaf*)leaf_array->At(6))->SetAddress(&nUnits);
        ((TLeaf*)leaf_array->At(7))->SetAddress(EdepTot);
      }
    };
    void WriteToASCII(const G4Event*);
    void WriteToROOT(const G4Event*);
    void readHitsCollection(const TrackerHitsCollection*, Det_event*);

    G4int SlowCollID, FastCollID, HighZCollID;
    G4int printModulo;
    
    Prim_event PrimEvent;
    Det_event SlowEvent, FastEvent, HighZEvent;
    
    std::ofstream* asciiFile;
    TFile* rootFile;
    TTree *G4dataTree;
};

#endif
