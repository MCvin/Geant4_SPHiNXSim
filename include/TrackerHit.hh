//
//  TrackerHit.hh 2013-09-04  Maxime Chauvin
//

#ifndef TrackerHit_h
#define TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"

class TrackerHit : public G4VHit
{
  public:
    TrackerHit();
   ~TrackerHit();
    TrackerHit(const TrackerHit&);
    const TrackerHit& operator=(const TrackerHit&);
    G4int operator==(const TrackerHit&) const;

    inline void *operator new(size_t);
    inline void operator delete(void* aHit);

  public:
    // set and get id of the scintillator
    void SetDaughterId(G4int id){ m_daughterid = id; };
    G4int GetDaughterId() { return m_daughterid; };
    void SetMotherId(G4int id){ m_motherid = id; };
    G4int GetMotherId() { return m_motherid; };

    // set and get energy deposition
    void SetEdep(G4double de) { m_Edep = de; };
    G4double GetEdep() { return m_Edep; };

    // set and get step length
    void SetStepLength(G4double sl) { m_stepLength = sl; };
    G4double GetStepLength() { return m_stepLength; };

    // set and get position
    void SetPosition(G4ThreeVector xyz) { m_position = xyz; };
    G4ThreeVector GetPosition() { return m_position; };

    // set and get post position
    void SetPostPosition(G4ThreeVector xyz) { m_postPosition = xyz; };
    G4ThreeVector GetPostPosition() { return m_postPosition; };

    // set and get momentum
    void SetMomentum(G4ThreeVector pxyz) { m_momentum = pxyz; };
    G4ThreeVector GetMomentum() { return m_momentum; };

    // set and get post momentum
    void SetPostMomentum(G4ThreeVector pxyz) { m_postMomentum = pxyz; };
    G4ThreeVector GetPostMomentum() { return m_postMomentum; };

    // set and get particle
    void SetParticle(G4DynamicParticle* par) { m_particle = par; };
    G4DynamicParticle* GetParticle() { return m_particle; };

    // set and get particle name/charge
    void SetParticleName(G4String name) { m_particleName = name; };
    G4String GetParticleName() { return m_particleName; };
    void SetCharge(G4double charge) { m_charge = charge; };
    G4double GetCharge() { return m_charge; };

    // set and get kinetic energy of the particle
    void SetKineticEnergy(G4double kine) { m_kineticEnergy = kine; };
    G4double GetKineticEnergy() { return m_kineticEnergy; };
    void SetPostKineticEnergy(G4double kine) { m_postKineticEnergy = kine; };
    G4double GetPostKineticEnergy() { return m_postKineticEnergy; };

    // set and get process name
    void SetProcessName(G4String name) { m_processName = name; };
    G4String GetProcessName() { return m_processName; };

    // set and get polarization vector
    void SetPolarization(G4ThreeVector pol) { m_polarization = pol; };
    G4ThreeVector GetPolarization() { return m_polarization; };
    void SetPostPolarization(G4ThreeVector ppol) { m_postPolarization = ppol; };
    G4ThreeVector GetPostPolarization() { return m_postPolarization; };

  private:
    // id of the scintillator
    G4int m_daughterid, m_motherid;
    // deposit energy and step length
    G4double m_Edep, m_stepLength; 
    // particle and particle name/charge
    G4DynamicParticle* m_particle;
    G4String m_particleName;
    G4double m_charge;
    // start and end position of the track
    G4ThreeVector m_position, m_postPosition;
    // start and end momentum of the track
    G4ThreeVector m_momentum, m_postMomentum;
    // kinetic energy of the particle
    G4double m_kineticEnergy, m_postKineticEnergy;
    // process name
    G4String m_processName;
    // polarization verctor
    G4ThreeVector m_polarization, m_postPolarization;
};

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;

extern G4Allocator<TrackerHit> TrackerHitAllocator;

inline void* TrackerHit::operator new(size_t)
{
  void* aHit = (void*)TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void TrackerHit::operator delete(void* aHit)
{
  TrackerHitAllocator.FreeSingle((TrackerHit*) aHit);
}

#endif
