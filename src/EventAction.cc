//
//  EventAction.cc 2013-09-04  Maxime Chauvin
//

#include "EventAction.hh"

#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(): SlowCollID(-1), FastCollID(-1), HighZCollID(-1)
{
  printModulo = 10000;
}

EventAction::~EventAction()
{
  asciiFile->close();
  G4dataTree->Write();
  rootFile->Close();
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  if (evt->GetEventID()%printModulo == 0){ 
    G4cout << "    ---> Begin of event: " << evt->GetEventID() << "\n" << G4endl;
  }

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if (SlowCollID==-1){SlowCollID = SDman->GetCollectionID("SlowSD");}
  if (FastCollID==-1){FastCollID = SDman->GetCollectionID("FastSD");}
  if (HighZCollID==-1){HighZCollID = SDman->GetCollectionID("HighZSD");}
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  // write to ROOT file
  WriteToROOT(evt);
  // write to ASCII file
  WriteToASCII(evt);
}

void EventAction::CreateROOTFile(const G4String name)
{
  rootFile = new TFile(name,"RECREATE");
  G4dataTree = new TTree( "G4data" , "PoGOSIM data" );
  G4dataTree->Branch("Primary", &PrimEvent, "Id/I:Energy/F:Px/F:Py/F:Pz/F:Dx/F:Dy/F:Dz/F");
  TBranch* Slow_branch = G4dataTree->Branch("Slow", &SlowEvent, "nhits/I:Edep[nhits]/F:Px[nhits]/F:Py[nhits]/F:Pz[nhits]/F:UnitId[nhits]/I:nUnits/I:EdepTot[60]/F");
  SlowEvent.SetAddress(Slow_branch); //Set the address for each leaf explicitly (IMPORTANT!)
  TBranch* Fast_branch = G4dataTree->Branch("Fast", &FastEvent, "nhits/I:Edep[nhits]/F:Px[nhits]/F:Py[nhits]/F:Pz[nhits]/F:UnitId[nhits]/I:nUnits/I:EdepTot[60]/F");
  FastEvent.SetAddress(Fast_branch); //Set the address for each leaf explicitly (IMPORTANT!)
  TBranch* HighZ_branch = G4dataTree->Branch("HighZ", &HighZEvent, "nhits/I:Edep[nhits]/F:Px[nhits]/F:Py[nhits]/F:Pz[nhits]/F:UnitId[nhits]/I:nUnits/I:EdepTot[60]/F");
  HighZEvent.SetAddress(HighZ_branch); //Set the address for each leaf explicitly (IMPORTANT!)
  G4dataTree->SetAutoSave(1000000000); // to avoid getting many trees in the root file
}

void EventAction::WriteToROOT(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  const G4PrimaryVertex* primaryVertex = evt->GetPrimaryVertex();
  G4PrimaryParticle* primary = primaryVertex->GetPrimary(0);

  // get everything we want to output
    // primary particule
    PrimEvent.Id = evt->GetEventID();
    PrimEvent.E = primary->GetKineticEnergy()/keV;
    PrimEvent.Px = primaryVertex->GetX0()/mm;
    PrimEvent.Py = primaryVertex->GetY0()/mm;
    PrimEvent.Pz = primaryVertex->GetZ0()/mm;
    PrimEvent.Dx = primary->GetPx();
    PrimEvent.Dy = primary->GetPy();
    PrimEvent.Dz = primary->GetPz();
    // Slow scintillator hits
    TrackerHitsCollection* SlowHC = 0;
    if (HCE){ SlowHC = (TrackerHitsCollection*)HCE->GetHC(SlowCollID); }
    if (SlowHC){ readHitsCollection(SlowHC, &SlowEvent); }
    // Fast scintillator hits
    TrackerHitsCollection* FastHC = 0;
    if (HCE){ FastHC = (TrackerHitsCollection*)HCE->GetHC(FastCollID); }
    if (FastHC){ readHitsCollection(FastHC, &FastEvent); }
    // HighZ scintillator hits
    TrackerHitsCollection* HighZHC = 0;
    if (HCE){ HighZHC = (TrackerHitsCollection*)HCE->GetHC(HighZCollID); }
    if (HighZHC){ readHitsCollection(HighZHC, &HighZEvent); }

  // fill the Tree only if something happened in the Slow, Fast or HighZ
  if ((SlowEvent.nhits > 0) || (FastEvent.nhits > 0) || (HighZEvent.nhits > 0)) G4dataTree->Fill();
}

void EventAction::readHitsCollection(const TrackerHitsCollection* HC, Det_event* detEvt)
{
  G4int nhits = 0;
  G4double Edep = 0;

  detEvt->nUnits = 0;
  for (G4int i = 0; i < 60; i++) detEvt->EdepTot[i] = 0.0;

  for (G4int i = 0; i < HC->entries(); i++){
    Edep = (*HC)[i]->GetEdep()/keV;
    if (Edep > 0){
      G4int UnitId = (*HC)[i]->GetMotherId() + (*HC)[i]->GetDaughterId();
      
      detEvt->Edep[nhits] = Edep;
      detEvt->Px[nhits] = (*HC)[i]->GetPostPosition().x()/mm;
      detEvt->Py[nhits] = (*HC)[i]->GetPostPosition().y()/mm;
      detEvt->Pz[nhits] = (*HC)[i]->GetPostPosition().z()/mm;
      detEvt->UnitId[nhits] = UnitId;
      nhits++;
      
      if (detEvt->EdepTot[UnitId] == 0) detEvt->nUnits += 1;
      detEvt->EdepTot[UnitId] += Edep;
    }
  }
  detEvt->nhits = nhits;

  return ;
}

void EventAction::WriteToASCII(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  G4int UnitId = 0;
  G4double TotEdep = 0;
  G4double* Edep = new G4double[88]; //number of detector units 4*7 + 2*30
  for (G4int i = 0; i < 88; i++) Edep[i] = 0.0;

  // print Fast scintillator
  TrackerHitsCollection* FastHC = 0;
  if (HCE){ FastHC = (TrackerHitsCollection*)HCE->GetHC(FastCollID); }
  for (G4int i = 0; i < FastHC->entries(); i++){
    UnitId = (*FastHC)[i]->GetMotherId() + (*FastHC)[i]->GetDaughterId();
    Edep[UnitId] += (*FastHC)[i]->GetEdep();
    TotEdep += (*FastHC)[i]->GetEdep();
  }

  // print HighZ scintillator
  TrackerHitsCollection* HighZHC = 0;
  if (HCE){ HighZHC = (TrackerHitsCollection*)HCE->GetHC(HighZCollID); }
    for (G4int i = 0; i < HighZHC->entries(); i++){
    UnitId = (*HighZHC)[i]->GetMotherId() + (*HighZHC)[i]->GetDaughterId();
    Edep[UnitId+28] += (*HighZHC)[i]->GetEdep();
    TotEdep += (*HighZHC)[i]->GetEdep();
  }
  
  if (TotEdep > 0){
    for (G4int i = 0; i < 88; i++){
      (*asciiFile) << Edep[i]/keV << " ";
    }
    (*asciiFile) << G4endl;
  }
  delete[] Edep;
}
