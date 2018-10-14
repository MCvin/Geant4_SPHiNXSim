//
//  DetectorConstruction.cc 2013-09-04  Maxime Chauvin
//

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4GDMLParser.hh"
#include "G4RunManager.hh"

DetectorConstruction::DetectorConstruction()
{
  // create commands for interactive definition of the detector
  detectorMessenger = new DetectorMessenger(this);

  // define the geometry parameters
  DefineGeometry();

  // define some rotation matrices
  DefineRotations();
}

DetectorConstruction::~DetectorConstruction()
{
  delete detectorMessenger;
}

void DetectorConstruction::DefineGeometry()
{
  SlowHeight = 0.0*mm;          //can be set

  FastOuter = 11.0*cm;          //can be set
  FastHeight = 6.0*cm;          //can be set
  FastCutHeight = 2.5*cm;       //can be set

  HighZThick = 0.4*cm;          //can be set
  HighZHeight = 5.9*cm;         //can be set (FastHeight+SlowHeight-CFRPThick_Ribs)

  ShieldHeight = 12.0*cm;       //can be set

  ShieldThick_Pb = 1.0*mm;      //can be set
  ShieldThick_Sn = 0.5*mm;      //can be set
  ShieldThick_Cu = 0.25*mm;     //can be set
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // to export the geometry
  G4GDMLParser parser;

  // Elements
  G4Element *elementH  = new G4Element("Hydrogen",     "H",   1.,    1.0079*g/mole);
  G4Element *elementC  = new G4Element("Carbon",       "C",   6.,    12.011*g/mole);
  G4Element *elementGd = new G4Element("Gadolinium",  "Gd",  64.,  157.2500*g/mole);
  G4Element *elementAl = new G4Element("Aluminium",   "Al",  13.,  26.98150*g/mole);
  G4Element *elementGa = new G4Element("Gallium",     "Ga",  31.,  69.72000*g/mole);
  G4Element *elementO  = new G4Element("Oxygen",       "O",   8.,  15.99940*g/mole);

  // materials
  G4NistManager * nist_manager = G4NistManager::Instance();
  G4Material * galactic = nist_manager->FindOrBuildMaterial("G4_Galactic");
  G4Material * Al = nist_manager->FindOrBuildMaterial("G4_Al");
  G4Material * Pb = nist_manager->FindOrBuildMaterial("G4_Pb");
  G4Material * Sn = nist_manager->FindOrBuildMaterial("G4_Sn");
  G4Material * Cu = nist_manager->FindOrBuildMaterial("G4_Cu");
  G4Material * Si = nist_manager->FindOrBuildMaterial("G4_Si");
  G4Material * plastic = nist_manager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  //G4Material * CsI = nist_manager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  //G4Material * NaI = nist_manager->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  G4Material * BGO = nist_manager->FindOrBuildMaterial("G4_BGO");
  G4Material* GaGG = new G4Material("GaGG", 6.63*g/cm3, 4);
  GaGG->AddElement(elementGd, 3);
  GaGG->AddElement(elementAl, 2);
  GaGG->AddElement(elementGa, 3);
  GaGG->AddElement(elementO, 12);

  // CRFP (Carbon Fiber Reinforced Polymer): M55 Quasiisotropic Layup
  G4Material * CFRP = new G4Material("CFRP", 6., 12.011*g/mole, 1.66*g/cm3);
  //Epoxy (for FR4)
  G4Material* Epoxy = new G4Material("Epoxy", 1.2*g/cm3, 2);
  Epoxy->AddElement(elementH, 2);
  Epoxy->AddElement(elementC, 2);
  //FR4 (Glass + Epoxy)
  G4Material * SiO2 = nist_manager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* FR4 = new G4Material("FR4", 1.86*g/cm3, 2);
  FR4->AddMaterial(SiO2, 0.528);
  FR4->AddMaterial(Epoxy, 0.472);


  // Geometry design
  const G4double Pi = std::acos(-1.0);
  G4double FastInner = FastOuter*std::cos(Pi/6);
  G4double FastFace = FastOuter/2;
  G4double FastPMWidth = FastOuter/4; //30.0*mm;
  G4double FastPMHeight = FastHeight/3; //22.0*mm;
  G4double CFRPThick_LowZ = 0.75*mm;
  G4double CFRPThick_HighZ = 0.5*mm;
  G4double CFRPThick_Ribs = 1.0*mm;
  G4double HighZPMWidth = 9.0*mm;
  G4double HighZPMLength = 10.6*mm;
  G4double HighZPMHeight = 1.65*mm;
  G4double ShieldHeight_CFRP = FastHeight+8.6*cm;
  G4double ShieldThick_CFRP = 3.0*mm;
  G4double ShieldTopThick_CFRP = 1.0*mm;
  G4double ShieldInner_Cu = FastInner*3+HighZThick*4+CFRPThick_HighZ*6+1.5*cm;
  G4double ShieldInner_Sn = ShieldInner_Cu+2*ShieldThick_Cu;
  G4double ShieldInner_Pb = ShieldInner_Sn+2*ShieldThick_Sn;
  G4double ShieldInner_CFRP = ShieldInner_Pb+2*ShieldThick_Pb;
  G4double ShieldOuter_CFRP = ShieldInner_CFRP+2*ShieldThick_CFRP;

  // WORLD
  const G4double WorldSize  = 1000.0*cm; 
  G4Box* s_World = new G4Box("World", WorldSize/2, WorldSize/2, WorldSize/2);
  G4LogicalVolume* l_World = new G4LogicalVolume(s_World, galactic, "World");
  G4VPhysicalVolume* p_World = new G4PVPlacement(0,	G4ThreeVector(), "World", l_World, NULL, false, 0);

  // SPHiNX mother volume
  //G4Box* s_Sphinx = new G4Box("Sphinx", 50.0*cm/2, 50.0*cm/2, 15.0*cm/2);
  G4Box* s_Sphinx = new G4Box("Sphinx", WorldSize/2, WorldSize/2, WorldSize/2);
  G4LogicalVolume* l_Sphinx = new G4LogicalVolume(s_Sphinx, galactic, "Sphinx");
  p_Sphinx = new G4PVPlacement(0, G4ThreeVector(), "Sphinx", l_Sphinx, p_World, false, 0);

  // Plastic detector mother volume
  const G4String PUnitName = "Plastic Unit";
  const G4double PUnitZ[2] = {0.0*cm,SlowHeight+FastHeight+FastPMHeight};
  const G4double PUnitInner[2] = {0.0*cm,0.0*cm};
  const G4double PUnitOuter[2] = {FastInner/2+CFRPThick_HighZ,FastInner/2+CFRPThick_HighZ};
  G4Polyhedra* s_PUnit = new G4Polyhedra(PUnitName, 30.*deg, 360.*deg, 6, 2, PUnitZ, PUnitInner, PUnitOuter);
  G4LogicalVolume* l_PUnit = new G4LogicalVolume(s_PUnit, galactic, PUnitName);

      G4Box* s_Box = new G4Box("s_Box", FastOuter/4, FastOuter/4, SlowHeight/2+FastHeight/2);

	  //// Slow (plastic)
	  //const G4double SHexZ[2] = {0.0*cm,SlowHeight};
	  //const G4double SHexInner[2] = {0.0*cm,0.0*cm};
	  //const G4double SHexOuter[2] = {FastInner/2,FastInner/2};
	  //G4Polyhedra* s_Slow = new G4Polyhedra("Slow", 30.*deg, 360.*deg, 6, 2, SHexZ, SHexInner, SHexOuter);
      //
      //G4IntersectionSolid* s1_Slow = new G4IntersectionSolid("Slow_1", s_Slow, s_Box, 0, G4ThreeVector(FastOuter/4+CFRPThick_LowZ/2,FastOuter/4+CFRPThick_LowZ/2,FastHeight/2));
      //G4IntersectionSolid* s2_Slow = new G4IntersectionSolid("Slow_2", s_Slow, s_Box, 0, G4ThreeVector(FastOuter/4+CFRPThick_LowZ/2,-(FastOuter/4+CFRPThick_LowZ/2),FastHeight/2));
      //G4IntersectionSolid* s3_Slow = new G4IntersectionSolid("Slow_3", s_Slow, s_Box, 0, G4ThreeVector(-(FastOuter/4+CFRPThick_LowZ/2),-(FastOuter/4+CFRPThick_LowZ/2),FastHeight/2));
      //G4IntersectionSolid* s4_Slow = new G4IntersectionSolid("Slow_4", s_Slow, s_Box, 0, G4ThreeVector(-(FastOuter/4+CFRPThick_LowZ/2),FastOuter/4+CFRPThick_LowZ/2,FastHeight/2));
	  //G4LogicalVolume* l1_Slow = new G4LogicalVolume(s1_Slow, plastic, "Slow_1");
	  //G4LogicalVolume* l2_Slow = new G4LogicalVolume(s2_Slow, plastic, "Slow_2");
	  //G4LogicalVolume* l3_Slow = new G4LogicalVolume(s3_Slow, plastic, "Slow_3");
	  //G4LogicalVolume* l4_Slow = new G4LogicalVolume(s4_Slow, plastic, "Slow_4");
	  //G4VPhysicalVolume* p1_Slow = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight+FastHeight), l1_Slow, "Slow_1", l_PUnit, false, 0);
	  //G4VPhysicalVolume* p2_Slow = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight+FastHeight), l2_Slow, "Slow_2", l_PUnit, false, 1);
	  //G4VPhysicalVolume* p3_Slow = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight+FastHeight), l3_Slow, "Slow_3", l_PUnit, false, 2);
	  //G4VPhysicalVolume* p4_Slow = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight+FastHeight), l4_Slow, "Slow_4", l_PUnit, false, 3);

	  // Fast (plastic)
	  const G4double FHexZ[2] = {0.0*cm,FastHeight};
	  const G4double FHexInner[2] = {0.0*cm,0.0*cm};
	  const G4double FHexOuter[2] = {FastInner/2,FastInner/2};
	  G4Polyhedra* s_FHex = new G4Polyhedra("s_Hex", 30.*deg, 360.*deg, 6, 2, FHexZ, FHexInner, FHexOuter);
	  const G4double FQuadZ[3] = {0.0*cm,FastCutHeight,FastHeight};
	  const G4double FQuadInner[3] = {0.0*cm,0.0*cm,0.0*cm};
	  const G4double FQuadOuter[3] = {FastPMWidth,FastOuter/2,FastOuter/2};
	  G4Polyhedra* s_FQuad = new G4Polyhedra("s_Quad", 45.*deg, 360.*deg, 4, 3, FQuadZ, FQuadInner, FQuadOuter);

	  G4IntersectionSolid* s_Fast = new G4IntersectionSolid("Fast", s_FHex, s_FQuad);
      G4IntersectionSolid* s1_Fast = new G4IntersectionSolid("Fast_1", s_Fast, s_Box, 0, G4ThreeVector(FastOuter/4+CFRPThick_LowZ/2,FastOuter/4+CFRPThick_LowZ/2,FastHeight/2));
      G4IntersectionSolid* s2_Fast = new G4IntersectionSolid("Fast_2", s_Fast, s_Box, 0, G4ThreeVector(FastOuter/4+CFRPThick_LowZ/2,-(FastOuter/4+CFRPThick_LowZ/2),FastHeight/2));
      G4IntersectionSolid* s3_Fast = new G4IntersectionSolid("Fast_3", s_Fast, s_Box, 0, G4ThreeVector(-(FastOuter/4+CFRPThick_LowZ/2),-(FastOuter/4+CFRPThick_LowZ/2),FastHeight/2));
      G4IntersectionSolid* s4_Fast = new G4IntersectionSolid("Fast_4", s_Fast, s_Box, 0, G4ThreeVector(-(FastOuter/4+CFRPThick_LowZ/2),FastOuter/4+CFRPThick_LowZ/2,FastHeight/2));
	  G4LogicalVolume* l1_Fast = new G4LogicalVolume(s1_Fast, plastic, "Fast_1");
	  G4LogicalVolume* l2_Fast = new G4LogicalVolume(s2_Fast, plastic, "Fast_2");
	  G4LogicalVolume* l3_Fast = new G4LogicalVolume(s3_Fast, plastic, "Fast_3");
	  G4LogicalVolume* l4_Fast = new G4LogicalVolume(s4_Fast, plastic, "Fast_4");
	  G4VPhysicalVolume* p1_Fast = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight), l1_Fast, "Fast_1", l_PUnit, false, 0);
	  G4VPhysicalVolume* p2_Fast = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight), l2_Fast, "Fast_2", l_PUnit, false, 1);
	  G4VPhysicalVolume* p3_Fast = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight), l3_Fast, "Fast_3", l_PUnit, false, 2);
	  G4VPhysicalVolume* p4_Fast = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,FastPMHeight), l4_Fast, "Fast_4", l_PUnit, false, 3);

	  // PMT
	  const G4String PMTName = "PMT";
	  G4Box* s_PMT = new G4Box(PMTName, FastPMWidth/2, FastPMWidth/2, FastPMHeight/2);
	  G4LogicalVolume* l_PMT = new G4LogicalVolume(s_PMT, Al, PMTName);
	  G4VPhysicalVolume* p1_PMT = new G4PVPlacement(0, G4ThreeVector(  FastPMWidth/2+CFRPThick_LowZ/2,   FastPMWidth/2+CFRPThick_LowZ/2, FastPMHeight/2), l_PMT, PMTName, l_PUnit, false, 0);
	  G4VPhysicalVolume* p2_PMT = new G4PVPlacement(0, G4ThreeVector(  FastPMWidth/2+CFRPThick_LowZ/2, -(FastPMWidth/2+CFRPThick_LowZ/2),FastPMHeight/2), l_PMT, PMTName, l_PUnit, false, 0);
	  G4VPhysicalVolume* p3_PMT = new G4PVPlacement(0, G4ThreeVector(-(FastPMWidth/2+CFRPThick_LowZ/2),-(FastPMWidth/2+CFRPThick_LowZ/2),FastPMHeight/2), l_PMT, PMTName, l_PUnit, false, 0);
	  G4VPhysicalVolume* p4_PMT = new G4PVPlacement(0, G4ThreeVector(-(FastPMWidth/2+CFRPThick_LowZ/2),  FastPMWidth/2+CFRPThick_LowZ/2, FastPMHeight/2), l_PMT, PMTName, l_PUnit, false, 0);

	  // CFRP structure
      G4Box* s_SheetX = new G4Box("s_SheetX", FastInner/2, CFRPThick_LowZ/2, SlowHeight/2+FastHeight/2+FastPMHeight/2);
      G4Box* s_SheetY = new G4Box("s_SheetY", CFRPThick_LowZ/2, FastOuter/2, SlowHeight/2+FastHeight/2+FastPMHeight/2);
      G4UnionSolid* s_CFRPSheet = new G4UnionSolid("CFRPSheet", s_SheetX, s_SheetY);
      const G4double CHexZ[2] = {0.0*cm,SlowHeight+FastHeight+FastPMHeight};
	  const G4double CHexInner[2] = {FastInner/2,FastInner/2};
	  const G4double CHexOuter[2] = {FastInner/2+CFRPThick_HighZ,FastInner/2+CFRPThick_HighZ};
	  G4Polyhedra* s_Honey = new G4Polyhedra("CFRPHoney", 30.*deg, 360.*deg, 6, 2, CHexZ, CHexInner, CHexOuter);
	  // Final shape
      const G4String CFRPLowZName = "CFRPLowZ";
	  G4UnionSolid* s_CFRPLowZ = new G4UnionSolid(CFRPLowZName, s_CFRPSheet, s_Honey, 0, G4ThreeVector(0.0*cm,0.0*cm,-(SlowHeight/2+FastHeight/2+FastPMHeight/2)));
	  G4LogicalVolume* l_CFRPLowZ = new G4LogicalVolume(s_CFRPLowZ, CFRP, CFRPLowZName);
	  G4VPhysicalVolume* p_CFRPLowZ = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,SlowHeight/2+FastHeight/2+FastPMHeight/2), l_CFRPLowZ, CFRPLowZName, l_PUnit, false, 0);

  // Place 7 PUnit
  const G4double x = FastInner + HighZThick + 2*CFRPThick_HighZ;
  const G4double y = x * std::sin(Pi/3);
  const G4double h = 15.0*cm/2 -SlowHeight -FastHeight -FastPMHeight -0.5*cm;
    // optionnal to get the position of each quadrant center
    G4double qx = (FastInner/2-CFRPThick_LowZ/2)/3*((FastOuter-CFRPThick_LowZ*3/2)/(FastFace*3/2-CFRPThick_LowZ));
    G4double qy = (FastFace/2-CFRPThick_LowZ/2)/4 + (FastFace-CFRPThick_LowZ/2)/4;
    //G4double qx = 21.49*mm; //23.21*mm for Plastic to highZ
    //G4double qy = 21.44*mm; //23.1*mm for Plastic to highZ
    PQ1 = G4ThreeVector( qx, qy,0.0*cm);
    PQ2 = G4ThreeVector( qx,-qy,0.0*cm);
    PQ3 = G4ThreeVector(-qx,-qy,0.0*cm);
    PQ4 = G4ThreeVector(-qx, qy,0.0*cm);
    PrintPositions(   0.*deg,    0,  0, h, "Plastic");
    PrintPositions(   0.*deg,    x,  0, h, "Plastic");
    PrintPositions(  60.*deg,  x/2, -y, h, "Plastic");
    PrintPositions( 120.*deg, -x/2, -y, h, "Plastic");
    PrintPositions( 180.*deg,   -x,  0, h, "Plastic");
    PrintPositions( 240.*deg, -x/2,  y, h, "Plastic");
    PrintPositions( 300.*deg,  x/2,  y, h, "Plastic");

  PlaceCopy(      0,    0,  0, h, l_PUnit, PUnitName, l_Sphinx,  0);
  PlaceCopy(      0,    x,  0, h, l_PUnit, PUnitName, l_Sphinx,  4);
  PlaceCopy( zRot60,  x/2, -y, h, l_PUnit, PUnitName, l_Sphinx,  8);
  PlaceCopy(zRot120, -x/2, -y, h, l_PUnit, PUnitName, l_Sphinx, 12);
  PlaceCopy(zRot180,   -x,  0, h, l_PUnit, PUnitName, l_Sphinx, 16);
  PlaceCopy(zRot240, -x/2,  y, h, l_PUnit, PUnitName, l_Sphinx, 20);
  PlaceCopy(zRot300,  x/2,  y, h, l_PUnit, PUnitName, l_Sphinx, 24);


  // HighZ detector mother volume
  const G4String CUnitName = "HighZ Unit";
  G4Box* s_CUnit = new G4Box(CUnitName, HighZThick/2, FastFace/2, HighZPMHeight/2+HighZHeight/2);
  G4LogicalVolume* l_CUnit = new G4LogicalVolume(s_CUnit, galactic, CUnitName);

  	  // HighZ
  	  const G4String HighZName = "HighZ";
  	  G4Box* s_HighZ = new G4Box(HighZName, HighZThick/2, FastFace/4-CFRPThick_LowZ/4, HighZHeight/2);
      G4LogicalVolume* l_HighZ = new G4LogicalVolume(s_HighZ, BGO, HighZName);
  	  G4VPhysicalVolume* p1_HighZ = new G4PVPlacement(0, G4ThreeVector(0.0*cm,  FastFace/4+CFRPThick_LowZ/4 ,HighZPMHeight/2), l_HighZ, HighZName, l_CUnit, false, 0);
  	  G4VPhysicalVolume* p2_HighZ = new G4PVPlacement(0, G4ThreeVector(0.0*cm,-(FastFace/4+CFRPThick_LowZ/4),HighZPMHeight/2), l_HighZ, HighZName, l_CUnit, false, 1);

	  // PM
	  const G4String PMName = "PM";
	  G4Box* s_PM = new G4Box(PMName, HighZPMWidth/2, HighZPMLength/2, HighZPMHeight/2);
	  G4LogicalVolume* l_PM = new G4LogicalVolume(s_PM, Al, PMName);
	  G4VPhysicalVolume* p1_PM = new G4PVPlacement(0, G4ThreeVector(0.0*cm,  FastFace/4+CFRPThick_LowZ/4 ,-HighZHeight/2), l_PM, PMName, l_CUnit, false, 0);
	  G4VPhysicalVolume* p2_PM = new G4PVPlacement(0, G4ThreeVector(0.0*cm,-(FastFace/4+CFRPThick_LowZ/4),-HighZHeight/2), l_PM, PMName, l_CUnit, false, 0);

  // Place 30 HighZ units
  const G4double x1 = FastInner/2 + HighZThick/2 + CFRPThick_HighZ;
  const G4double y1 = x1 * std::sin(Pi/3);
  const G4double h1 = 15.0*cm/2 - HighZPMHeight/2 - HighZHeight/2 -CFRPThick_Ribs - 0.5*cm;
    // optionnal to get the position of each highZ center
    HZ1 = G4ThreeVector( 0.0*cm,  FastFace/4+CFRPThick_LowZ/4 ,0.0*cm);
    HZ2 = G4ThreeVector( 0.0*cm,-(FastFace/4+CFRPThick_LowZ/4),0.0*cm);
    PrintPositions(   0.*deg,          x1,   0, h1, "HighZ");
    PrintPositions(  60.*deg,        x1/2, -y1, h1, "HighZ");
    PrintPositions( 120.*deg,       -x1/2, -y1, h1, "HighZ");
    PrintPositions( 180.*deg,         -x1,   0, h1, "HighZ");
    PrintPositions( 240.*deg,       -x1/2,  y1, h1, "HighZ");
    PrintPositions( 300.*deg,        x1/2,  y1, h1, "HighZ");
    PrintPositions( 240.*deg,    -x1/2+x,   y1, h1, "HighZ");
    PrintPositions( 300.*deg,    -x1/2+x,  -y1, h1, "HighZ");
    PrintPositions(   0.*deg,    -x1+x/2,  0-y, h1, "HighZ");
    PrintPositions(  60.*deg,  -x1/2-x/2, y1-y, h1, "HighZ");
    PrintPositions( 120.*deg,     x1/2-x,   y1, h1, "HighZ");
    PrintPositions( 180.*deg,     x1-x/2,  0+y, h1, "HighZ");    
    PrintPositions(   0.*deg,      x1+x,     0, h1, "HighZ");
    PrintPositions(  60.*deg,    x1/2+x,   -y1, h1, "HighZ");
    PrintPositions(   0.*deg,    x1+x/2,   0-y, h1, "HighZ");
    PrintPositions(  60.*deg,  x1/2+x/2, -y1-y, h1, "HighZ");
    PrintPositions( 120.*deg, -x1/2+x/2, -y1-y, h1, "HighZ");
    PrintPositions(  60.*deg,  x1/2-x/2, -y1-y, h1, "HighZ");
    PrintPositions( 120.*deg, -x1/2-x/2, -y1-y, h1, "HighZ");
    PrintPositions( 180.*deg,   -x1-x/2,   0-y, h1, "HighZ");
    PrintPositions( 120.*deg,   -x1/2-x,   -y1, h1, "HighZ");
    PrintPositions( 180.*deg,     -x1-x,     0, h1, "HighZ");
    PrintPositions( 240.*deg,   -x1/2-x,    y1, h1, "HighZ");
    PrintPositions( 180.*deg,   -x1-x/2,   0+y, h1, "HighZ");
    PrintPositions( 240.*deg, -x1/2-x/2,  y1+y, h1, "HighZ");
    PrintPositions( 300.*deg,  x1/2-x/2,  y1+y, h1, "HighZ");
    PrintPositions( 240.*deg, -x1/2+x/2,  y1+y, h1, "HighZ");
    PrintPositions( 300.*deg,  x1/2+x/2,  y1+y, h1, "HighZ");
    PrintPositions(   0.*deg,    x1+x/2,   0+y, h1, "HighZ");
    PrintPositions( 300.*deg,    x1/2+x,    y1, h1, "HighZ");

    // ring 1
    PlaceCopy(      0,    x1,   0, h1, l_CUnit, CUnitName, l_Sphinx,  0);
    PlaceCopy( zRot60,  x1/2, -y1, h1, l_CUnit, CUnitName, l_Sphinx,  2);
    PlaceCopy(zRot120, -x1/2, -y1, h1, l_CUnit, CUnitName, l_Sphinx,  4);
    PlaceCopy(zRot180,   -x1,   0, h1, l_CUnit, CUnitName, l_Sphinx,  6);
    PlaceCopy(zRot240, -x1/2,  y1, h1, l_CUnit, CUnitName, l_Sphinx,  8);
    PlaceCopy(zRot300,  x1/2,  y1, h1, l_CUnit, CUnitName, l_Sphinx, 10);
    // ring 2
    PlaceCopy(zRot240,   -x1/2+x,   y1, h1, l_CUnit, CUnitName, l_Sphinx, 12);
    PlaceCopy(zRot300,   -x1/2+x,  -y1, h1, l_CUnit, CUnitName, l_Sphinx, 14);
    PlaceCopy(      0,   -x1+x/2,  0-y, h1, l_CUnit, CUnitName, l_Sphinx, 16);
    PlaceCopy( zRot60, -x1/2-x/2, y1-y, h1, l_CUnit, CUnitName, l_Sphinx, 18);
    PlaceCopy(zRot120,    x1/2-x,   y1, h1, l_CUnit, CUnitName, l_Sphinx, 20);
    PlaceCopy(zRot180,    x1-x/2,  0+y, h1, l_CUnit, CUnitName, l_Sphinx, 22);    
    // ring 3
    PlaceCopy(      0,      x1+x,     0, h1, l_CUnit, CUnitName, l_Sphinx, 24);
    PlaceCopy( zRot60,    x1/2+x,   -y1, h1, l_CUnit, CUnitName, l_Sphinx, 26);
    PlaceCopy(      0,    x1+x/2,   0-y, h1, l_CUnit, CUnitName, l_Sphinx, 28);
    PlaceCopy( zRot60,  x1/2+x/2, -y1-y, h1, l_CUnit, CUnitName, l_Sphinx, 30);
    PlaceCopy(zRot120, -x1/2+x/2, -y1-y, h1, l_CUnit, CUnitName, l_Sphinx, 32);
    PlaceCopy( zRot60,  x1/2-x/2, -y1-y, h1, l_CUnit, CUnitName, l_Sphinx, 34);
    PlaceCopy(zRot120, -x1/2-x/2, -y1-y, h1, l_CUnit, CUnitName, l_Sphinx, 36);
    PlaceCopy(zRot180,   -x1-x/2,   0-y, h1, l_CUnit, CUnitName, l_Sphinx, 38);
    PlaceCopy(zRot120,   -x1/2-x,   -y1, h1, l_CUnit, CUnitName, l_Sphinx, 40);
    PlaceCopy(zRot180,     -x1-x,     0, h1, l_CUnit, CUnitName, l_Sphinx, 42);
    PlaceCopy(zRot240,   -x1/2-x,    y1, h1, l_CUnit, CUnitName, l_Sphinx, 44);
    PlaceCopy(zRot180,   -x1-x/2,   0+y, h1, l_CUnit, CUnitName, l_Sphinx, 46);
    PlaceCopy(zRot240, -x1/2-x/2,  y1+y, h1, l_CUnit, CUnitName, l_Sphinx, 48);
    PlaceCopy(zRot300,  x1/2-x/2,  y1+y, h1, l_CUnit, CUnitName, l_Sphinx, 50);
    PlaceCopy(zRot240, -x1/2+x/2,  y1+y, h1, l_CUnit, CUnitName, l_Sphinx, 52);
    PlaceCopy(zRot300,  x1/2+x/2,  y1+y, h1, l_CUnit, CUnitName, l_Sphinx, 54);
    PlaceCopy(      0,    x1+x/2,   0+y, h1, l_CUnit, CUnitName, l_Sphinx, 56);
    PlaceCopy(zRot300,    x1/2+x,    y1, h1, l_CUnit, CUnitName, l_Sphinx, 58);

  // some passive volumes to represent the electronics
    // analog board
    const G4String AnaName = "Analog";
    G4Tubs* s_Ana = new G4Tubs(AnaName, 0.0*cm, ShieldInner_Cu/2, 3.0*mm/2, 0.*deg, 360.*deg);
    G4LogicalVolume* l_Ana = new G4LogicalVolume(s_Ana, FR4, AnaName);
    G4VPhysicalVolume* p_Ana = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,h-3.0*mm/2), l_Ana, AnaName, l_Sphinx, false, 0);
    // digital board
    const G4String DigiName = "Digital";
    G4Tubs* s_Digi = new G4Tubs(DigiName, 0.0*cm, ShieldInner_Cu/2, 2.0*mm/2, 0.*deg, 360.*deg);
    G4LogicalVolume* l_Digi = new G4LogicalVolume(s_Digi, FR4, DigiName);
    G4VPhysicalVolume* p_Digi = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,h-2.5*cm), l_Digi, DigiName, l_Sphinx, false, 0);
    //// HV PSU
    //const G4String HVName = "HV PSU";
    //G4Box* s_HV = new G4Box(HVName, 10.0*cm/2, 10.0*cm/2, 8.0*cm/2);
    //G4LogicalVolume* l_HV = new G4LogicalVolume(s_HV, Al, HVName);
    //G4VPhysicalVolume* p_HV = new G4PVPlacement(0, G4ThreeVector(19.0*cm,19.0*cm,-2.5*cm), l_HV, HVName, l_Sphinx, false, 0);
    //// LV PSU
    //const G4String LVName = "LV PSU";
    //G4Box* s_LV = new G4Box(LVName, 10.0*cm/2, 10.0*cm/2, 8.0*cm/2);
    //G4LogicalVolume* l_LV = new G4LogicalVolume(s_LV, Al, LVName);
    //G4VPhysicalVolume* p_LV = new G4PVPlacement(0, G4ThreeVector(-19.0*cm,-19.0*cm,-2.5*cm), l_LV, LVName, l_Sphinx, false, 0);
    //// HV1 PSU
    //const G4String HV1Name = "HV1 PSU";
    //G4Box* s_HV1 = new G4Box(HV1Name, 10.0*cm/2, 10.0*cm/2, 8.0*cm/2);
    //G4LogicalVolume* l_HV1 = new G4LogicalVolume(s_HV1, Al, HV1Name);
    //G4VPhysicalVolume* p_HV1 = new G4PVPlacement(0, G4ThreeVector(19.0*cm,-19.0*cm,-2.5*cm), l_HV1, HV1Name, l_Sphinx, false, 0);
    //// LV1 PSU
    //const G4String LV1Name = "LV1 PSU";
    //G4Box* s_LV1 = new G4Box(LV1Name, 10.0*cm/2, 10.0*cm/2, 8.0*cm/2);
    //G4LogicalVolume* l_LV1 = new G4LogicalVolume(s_LV1, Al, LV1Name);
    //G4VPhysicalVolume* p_LV1 = new G4PVPlacement(0, G4ThreeVector(-19.0*cm,19.0*cm,-2.5*cm), l_LV1, LV1Name, l_Sphinx, false, 0);

  // Shield Pb + Sn + Cu (photon shield)
    const G4String ShieldName_Pb = "Shield_Pb";
    const G4String ShieldName_Sn = "Shield_Sn";
    const G4String ShieldName_Cu = "Shield_Cu";
    G4Tubs* s_Shield_Pb = new G4Tubs(ShieldName_Pb, ShieldInner_Pb/2, ShieldInner_Pb/2+ShieldThick_Pb, ShieldHeight/2, 0.*deg, 360.*deg);
    G4Tubs* s_Shield_Sn = new G4Tubs(ShieldName_Sn, ShieldInner_Sn/2, ShieldInner_Sn/2+ShieldThick_Sn, ShieldHeight/2, 0.*deg, 360.*deg);
    G4Tubs* s_Shield_Cu = new G4Tubs(ShieldName_Cu, ShieldInner_Cu/2, ShieldInner_Cu/2+ShieldThick_Cu, ShieldHeight/2, 0.*deg, 360.*deg);
    G4LogicalVolume* l_Shield_Pb = new G4LogicalVolume(s_Shield_Pb, Pb, ShieldName_Pb);
    G4LogicalVolume* l_Shield_Sn = new G4LogicalVolume(s_Shield_Sn, Sn, ShieldName_Sn);
    G4LogicalVolume* l_Shield_Cu = new G4LogicalVolume(s_Shield_Cu, Cu, ShieldName_Cu);
    G4VPhysicalVolume* p_Shield_Pb = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,15.0*cm/2-ShieldHeight/2-0.5*cm), l_Shield_Pb, ShieldName_Pb, l_Sphinx, false, 0);
    G4VPhysicalVolume* p_Shield_Sn = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,15.0*cm/2-ShieldHeight/2-0.5*cm), l_Shield_Sn, ShieldName_Sn, l_Sphinx, false, 0);
    G4VPhysicalVolume* p_Shield_Cu = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,15.0*cm/2-ShieldHeight/2-0.5*cm), l_Shield_Cu, ShieldName_Cu, l_Sphinx, false, 0);

  // Shield CFRP (charged particle shield)
    const G4String ShieldName_CFRP = "Shield_CFRP";
    const G4double Shield_CFRPZ[4] = {0.*cm,ShieldHeight_CFRP-ShieldTopThick_CFRP,ShieldHeight_CFRP-ShieldTopThick_CFRP,ShieldHeight_CFRP};
    const G4double Shield_CFRPInner[4] = {ShieldInner_CFRP/2,ShieldInner_CFRP/2,0.0*cm,0.0*cm};
    const G4double Shield_CFRPOuter[4] = {ShieldOuter_CFRP/2,ShieldOuter_CFRP/2,ShieldOuter_CFRP/2,ShieldOuter_CFRP/2};
    G4Polycone* s_Cyl_CFRP = new G4Polycone(ShieldName_CFRP, 0.*deg, 360.*deg, 4, Shield_CFRPZ, Shield_CFRPInner, Shield_CFRPOuter);
    // Shield ribs
	const G4double RHexZ[2] = {-0.1*cm,CFRPThick_Ribs};
	const G4double RHexInner[2] = {0.0*cm,0.0*cm};
	const G4double RHexOuter[2] = {FastInner/2,FastInner/2};
	G4Polyhedra* s_CFRPRib1 = new G4Polyhedra(ShieldName_CFRP, 30.*deg, 360.*deg, 6, 2, RHexZ, RHexInner, RHexOuter);
    G4UnionSolid* s_CFRPRib2 = new G4UnionSolid(ShieldName_CFRP, s_CFRPRib1, s_CFRPRib1, 0, G4ThreeVector(-x,0.0*cm,0.0*cm));
    G4UnionSolid* s_CFRPRib4 = new G4UnionSolid(ShieldName_CFRP, s_CFRPRib2, s_CFRPRib2, 0, G4ThreeVector(x/2,-y,0.0*cm));
    G4UnionSolid* s_CFRPRib6 = new G4UnionSolid(ShieldName_CFRP, s_CFRPRib4, s_CFRPRib2, 0, G4ThreeVector(x/2,y,0.0*cm));
    G4UnionSolid* s_CFRPRib7 = new G4UnionSolid(ShieldName_CFRP, s_CFRPRib6, s_CFRPRib1, 0, G4ThreeVector(x,0.0*cm,0.0*cm));
    G4Tubs* s_CFRPRibCyl = new G4Tubs(ShieldName_CFRP, 0.0*cm, ShieldInner_CFRP/2-(ShieldThick_Pb+ShieldThick_Sn+ShieldThick_Cu), CFRPThick_Ribs/2, 0.*deg, 360.*deg);
    G4SubtractionSolid* s_CFRPRibs = new G4SubtractionSolid(ShieldName_CFRP, s_CFRPRibCyl, s_CFRPRib7, 0, G4ThreeVector(0.0*cm,0.0*cm,-CFRPThick_Ribs/2));
    // Final shape
    G4UnionSolid* s_Shield_CFRP = new G4UnionSolid(ShieldName_CFRP, s_Cyl_CFRP, s_CFRPRibs, 0, G4ThreeVector(0.0*cm,0.0*cm,ShieldHeight_CFRP-ShieldTopThick_CFRP-CFRPThick_Ribs/2));
    G4LogicalVolume* l_Shield_CFRP = new G4LogicalVolume(s_Shield_CFRP, CFRP, ShieldName_CFRP);
    G4VPhysicalVolume* p_Shield_CFRP = new G4PVPlacement(0, G4ThreeVector(0.0*cm,0.0*cm,15.0*cm/2-ShieldHeight_CFRP+ShieldTopThick_CFRP-0.5*cm), l_Shield_CFRP, ShieldName_CFRP, l_Sphinx, false, 0);

  // Export geometry in GDML format.
    //parser.Write("SPHiNXnew.gdml", p_Sphinx,false);

  // sensitive detectors
  G4SDManager *SDman = G4SDManager::GetSDMpointer();

    SlowSD = (TrackerSD*)SDman->FindSensitiveDetector("SlowSD");
    if (SlowSD) delete SlowSD; //G4SDManager::Activate("SlowSD",false)
      SlowSD = new TrackerSD("SlowSD");
      SDman->AddNewDetector(SlowSD);
      //l1_Slow->SetSensitiveDetector(SlowSD);
      //l2_Slow->SetSensitiveDetector(SlowSD);
      //l3_Slow->SetSensitiveDetector(SlowSD);
      //l4_Slow->SetSensitiveDetector(SlowSD);

    FastSD = (TrackerSD*)SDman->FindSensitiveDetector("FastSD");
    if (FastSD) delete FastSD;
      FastSD = new TrackerSD("FastSD");
      SDman->AddNewDetector(FastSD);
      l1_Fast->SetSensitiveDetector(FastSD);
      l2_Fast->SetSensitiveDetector(FastSD);
      l3_Fast->SetSensitiveDetector(FastSD);
      l4_Fast->SetSensitiveDetector(FastSD);

    HighZSD = (TrackerSD*)SDman->FindSensitiveDetector("HighZSD");
    if (HighZSD) delete HighZSD;
      HighZSD = new TrackerSD("HighZSD");
      SDman->AddNewDetector(HighZSD);
      l_HighZ->SetSensitiveDetector(HighZSD);


  // visualization
  G4VisAttributes* colourBlack = new G4VisAttributes(G4Colour(0.,0.,0.));
  G4VisAttributes* colourWhite = new G4VisAttributes(G4Colour(1.,1.,1.));
  G4VisAttributes* colourTWhite = new G4VisAttributes(G4Colour(1.,1.,1.,0.2));
  G4VisAttributes* colourBlue = new G4VisAttributes(G4Colour(0.,0.,1.,0.4));
  G4VisAttributes* colourLightBlue = new G4VisAttributes(G4Colour(0.5,0.6,1.,0.4));
  G4VisAttributes* colourGreen = new G4VisAttributes(G4Colour(0.,1.,0.));
  G4VisAttributes* colourGrey = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  G4VisAttributes* colourDarkGrey = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  G4VisAttributes* colourRed = new G4VisAttributes(G4Colour(1.,0.,0.,0.4));
  colourTWhite->SetForceSolid(true);
  colourWhite->SetForceSolid(true);
  colourBlue->SetForceSolid(true);
  colourLightBlue->SetForceSolid(true);
  //colourGreen->SetForceSolid(true);
  colourGrey->SetForceSolid(true);
  colourDarkGrey->SetForceSolid(true);
  colourRed->SetForceSolid(true);

  //l1_Slow->SetVisAttributes(colourLightBlue);
  //l2_Slow->SetVisAttributes(colourLightBlue);
  //l3_Slow->SetVisAttributes(colourLightBlue);
  //l4_Slow->SetVisAttributes(colourLightBlue);
  l1_Fast->SetVisAttributes(colourBlue);
  l2_Fast->SetVisAttributes(colourBlue);
  l3_Fast->SetVisAttributes(colourBlue);
  l4_Fast->SetVisAttributes(colourBlue);
  l_CFRPLowZ->SetVisAttributes(colourWhite);
  l_HighZ->SetVisAttributes(colourRed);
  l_PMT->SetVisAttributes(colourGrey);
  l_PM->SetVisAttributes(colourDarkGrey);
  l_Ana->SetVisAttributes(colourGreen);
  l_Digi->SetVisAttributes(colourGreen);
  //l_HV->SetVisAttributes(colourBlack);
  //l_HV1->SetVisAttributes(colourBlack);
  //l_LV->SetVisAttributes(colourBlack);
  //l_LV1->SetVisAttributes(colourBlack);
  //l_Shield_Pb->SetVisAttributes(colourGrey);
  //l_Shield_Sn->SetVisAttributes(colourLightBlue);
  //l_Shield_Cu->SetVisAttributes(colourBlue);
  l_Shield_CFRP->SetVisAttributes(colourTWhite);
  l_Shield_Pb->SetVisAttributes(colourDarkGrey);
  l_Shield_Sn->SetVisAttributes(colourGrey);
  l_Shield_Cu->SetVisAttributes(colourLightBlue);

  l_World->SetVisAttributes(G4VisAttributes::Invisible);
  l_Sphinx->SetVisAttributes(G4VisAttributes::Invisible);
  //l_Sphinx->SetVisAttributes(colourBlack);
  l_PUnit->SetVisAttributes(G4VisAttributes::Invisible);
  l_CUnit->SetVisAttributes(G4VisAttributes::Invisible);

  return p_World;
}

void DetectorConstruction::PlaceCopy(G4RotationMatrix* R, G4double X, G4double Y, G4double Z, G4LogicalVolume* Child, 
									 G4String Name, G4LogicalVolume* Parent, G4int CopyNum)
{
  new G4PVPlacement(R, G4ThreeVector(X*mm,Y*mm,Z*mm), Child, Name, Parent, false, CopyNum);
  return ;
}

void DetectorConstruction::PrintPositions(G4double rotZ, G4double X, G4double Y, G4double Z, G4String Name)
{
  if (Name == "Plastic") {
    G4ThreeVector P1 = PQ1;
    G4ThreeVector P2 = PQ2;
    G4ThreeVector P3 = PQ3;
    G4ThreeVector P4 = PQ4;
    P1 = P1.rotateZ(-rotZ);
    P2 = P2.rotateZ(-rotZ);
    P3 = P3.rotateZ(-rotZ);
    P4 = P4.rotateZ(-rotZ);
    G4cout << X+P1.getX() <<" "<< Y+P1.getY() << G4endl;
    G4cout << X+P2.getX() <<" "<< Y+P2.getY() << G4endl;
    G4cout << X+P3.getX() <<" "<< Y+P3.getY() << G4endl;
    G4cout << X+P4.getX() <<" "<< Y+P4.getY() << G4endl;
  }
  if (Name == "HighZ") {
    G4ThreeVector P1 = HZ1;
    G4ThreeVector P2 = HZ2;
    P1 = P1.rotateZ(-rotZ);
    P2 = P2.rotateZ(-rotZ);
    G4cout << X+P1.getX() <<" "<< Y+P1.getY() << G4endl;
    G4cout << X+P2.getX() <<" "<< Y+P2.getY() << G4endl;
  }

  return ;
}

void DetectorConstruction::SetInstRotY(G4double alpha)
{
  InstRotY = alpha;
  *yRotInst=*Rot0;
  yRotInst->rotateY(InstRotY);
  // rotate the mother volume Sphinx to rotate all the instrument
  p_Sphinx->SetRotation(yRotInst);

  // update the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "  Rotation Y of the polarimeter set to " << alpha/deg << " deg" << G4endl;
}
void DetectorConstruction::SetFastOuter(G4double length)
{
  FastOuter = length;
  UpdateGeometry();
  G4cout << "  Plastic size set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::SetFastHeight(G4double length)
{
  FastHeight = length;
  UpdateGeometry();
  G4cout << "  Plastic height set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::SetHighZThick(G4double length)
{
  HighZThick = length;
  UpdateGeometry();
  G4cout << "  HighZ thickness set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::SetHighZHeight(G4double length)
{
  HighZHeight = length;
  UpdateGeometry();
  G4cout << "  HighZ height set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::SetShieldHeight(G4double length)
{
  ShieldHeight = length;
  UpdateGeometry();
  G4cout << "  Shield height set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::SetShieldThick_Pb(G4double length)
{
  ShieldThick_Pb = length;
  UpdateGeometry();
  G4cout << "  Shield Pb thickness set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::SetShieldThick_Sn(G4double length)
{
  ShieldThick_Sn = length;
  UpdateGeometry();
  G4cout << "  Shield Sn thickness set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::SetShieldThick_Cu(G4double length)
{
  ShieldThick_Cu = length;
  UpdateGeometry();
  G4cout << "  Shield Cu thickness set to " << length/cm << " cm" << G4endl;
}
void DetectorConstruction::UpdateGeometry()
{
  // clean-up previous geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  //define new one
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4cout << "Geometry has been updated!" << G4endl;
}
