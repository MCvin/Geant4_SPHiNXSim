#ifndef G4data_h
#define G4data_h
#include <TTree.h>
#include <TLeaf.h>

struct Prim_event{
  int Id;
  float E;
  float Px, Py, Pz;
  float Dx, Dy, Dz;
  inline void SetAddress(TBranch* branch){
    branch->SetAddress(this);
  }
};
struct Det_event{
  int nhits;
  float Edep[1000];
  float Px[1000];
  float Py[1000];
  float Pz[1000];
  int UnitId[1000];
  int nUnits;
  float EdepTot[60];
  inline void SetAddress(TBranch* branch){
    branch->SetAddress(this);
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
struct G4data{
  Prim_event* Primary;
  Det_event* Slow;
  Det_event* Fast;
  Det_event* HighZ;
  G4data(){
  Primary = new Prim_event;
  Slow = new Det_event;
  Fast = new Det_event;
  HighZ = new Det_event;
  }
  inline void SetAddress(TTree* tree){
    Primary->SetAddress(tree->GetBranch("Primary"));
    Slow->SetAddress(tree->GetBranch("Slow"));
    Fast->SetAddress(tree->GetBranch("Fast"));
    HighZ->SetAddress(tree->GetBranch("HighZ"));
  }
};

#endif
