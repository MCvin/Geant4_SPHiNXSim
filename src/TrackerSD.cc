//
//  TrackerSD.cc 2013-09-04  Maxime Chauvin
//

#include "TrackerSD.hh"

#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"

TrackerSD::TrackerSD(G4String name) : G4VSensitiveDetector(name)
{
  collectionName.insert(name);
}

TrackerSD::~TrackerSD(){}

void TrackerSD::Initialize(G4HCofThisEvent* )
{
  // create a new "hits collection"
  m_HitsCollection = new TrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
}

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // get track, deposit energy and step length
  G4Track* aTrack = aStep->GetTrack();
  G4double Edep = aStep->GetTotalEnergyDeposit();
  G4double stepLength = aStep->GetStepLength();

  // if the process is tranportation, do nothing
  G4String processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  if (processName=="Transportation" && Edep==0.) {
    return true;
  }

  // get scintillator id
  G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4int Daughterid = theTouchable->GetCopyNumber(0); //CopyNumber of the detector daughter volume
  G4int Motherid = theTouchable->GetCopyNumber(1);	 //CopyNumber of the detector mother volume

  // get particle name
  G4String particleName = aTrack->GetDefinition()->GetParticleName();

  // create a new "hit"
  TrackerHit* newHit = new TrackerHit();

  // set energy deposition, step length, current and post position,
  // particle name, total energy of the particle, process name
  // and current/post polarization vector.
  newHit->SetCharge( aTrack->GetDefinition()->GetPDGCharge());
  newHit->SetDaughterId( Daughterid );
  newHit->SetMotherId( Motherid );
  newHit->SetEdep( Edep );
  newHit->SetStepLength( stepLength );
  newHit->SetPosition( aStep->GetPreStepPoint()->GetPosition() );
  newHit->SetPostPosition( aStep->GetPostStepPoint()->GetPosition() );
  newHit->SetMomentum( aStep->GetPreStepPoint()->GetMomentumDirection() );
  newHit->SetPostMomentum( aStep->GetPostStepPoint()->GetMomentumDirection() );
  newHit->SetParticleName( particleName );
  newHit->SetKineticEnergy( aStep->GetPreStepPoint()->GetKineticEnergy() );
  newHit->SetPostKineticEnergy( aStep->GetPostStepPoint()->GetKineticEnergy() );
  newHit->SetProcessName( aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
  newHit->SetPolarization( aStep->GetPreStepPoint()->GetPolarization() );
  newHit->SetPostPolarization( aStep->GetPostStepPoint()->GetPolarization() );
  
  // register "hit"(newHit) to "hits collection"
  m_HitsCollection->insert(newHit);

  return true;
}

void TrackerSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  // register "hits collection" to HCE
  HCE->AddHitsCollection(HCID, m_HitsCollection); 
}
