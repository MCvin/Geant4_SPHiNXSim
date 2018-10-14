//
//  DetectorConstruction.hh 2013-09-04  Maxime Chauvin
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
   ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

    void SetInstRotY(G4double);
    inline G4double GetInstRotY() { return InstRotY; };
    void SetFastOuter(G4double);
    void SetFastHeight(G4double);
    void SetHighZThick(G4double);
    void SetHighZHeight(G4double);
    void SetShieldHeight(G4double);
    void SetShieldThick_Pb(G4double);
    void SetShieldThick_Sn(G4double);
    void SetShieldThick_Cu(G4double);

  private:
    DetectorMessenger* detectorMessenger;
    G4VPhysicalVolume* p_Sphinx;
    
    G4double SlowHeight;
    G4double FastOuter;
    G4double FastHeight;
    G4double FastCutHeight;
    G4double HighZThick;
    G4double HighZHeight;
    G4double ShieldHeight;
    G4double ShieldThick_Pb;
    G4double ShieldThick_Sn;
    G4double ShieldThick_Cu;
    G4double InstRotY;
    //Some rotation matrices
    G4RotationMatrix* Rot0;
    G4RotationMatrix* yRotInst;
    G4RotationMatrix* zRot60;
    G4RotationMatrix* zRot120;
    G4RotationMatrix* zRot180;
    G4RotationMatrix* zRot240;
    G4RotationMatrix* zRot300;

    G4ThreeVector PQ1, PQ2, PQ3, PQ4;
    G4ThreeVector HZ1, HZ2;

    TrackerSD* SlowSD;
    TrackerSD* FastSD;
    TrackerSD* HighZSD;

    void DefineGeometry();
    void DefineRotations()
    {
      Rot0 = new G4RotationMatrix();
      yRotInst = new G4RotationMatrix();
      zRot60 = new G4RotationMatrix();
      zRot60->rotateZ(60.*deg);
      zRot120 = new G4RotationMatrix();
      zRot120->rotateZ(120.*deg);
      zRot180 = new G4RotationMatrix();
      zRot180->rotateZ(180.*deg);
      zRot240 = new G4RotationMatrix();
      zRot240->rotateZ(240.*deg);
      zRot300 = new G4RotationMatrix();
      zRot300->rotateZ(300.*deg);
	}
	void PlaceCopy(G4RotationMatrix*, G4double, G4double, G4double, G4LogicalVolume*, 
				   G4String, G4LogicalVolume*, G4int);
    void PrintPositions(G4double, G4double, G4double, G4double, G4String);
    void UpdateGeometry();

};

#endif
