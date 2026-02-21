#ifndef RIPTIDE_PHOTOCATHODE_HIT_HPP
#define RIPTIDE_PHOTOCATHODE_HIT_HPP

#include <G4Allocator.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

namespace riptide
{

  // Classe PhotocathodeHit deriva da G4VHit e rappresenta un singolo hit su un fotocatodo.
  class PhotocathodeHit : public G4VHit
  {
    G4ThreeVector m_global_position; // Posizione del fotocatodo colpito (mm)
    double m_time{};                 // Tempo globale (ns)
    double m_energy{};               // Energia del fotone (eV)
    int m_detector_id{};             // 0 = X, 1 = Y, etc.

  public:
    PhotocathodeHit() = default;

    void set_global_position(G4ThreeVector const &pos)
    {
      m_global_position = pos;
    }
    void set_time(double t)
    {
      m_time = t;
    }
    void set_energy(double e)
    {
      m_energy = e;
    }
    void set_detector_id(int id)
    {
      m_detector_id = id;
    }

    G4ThreeVector const &global_position() const
    {
      return m_global_position;
    }
    double time() const
    {
      return m_time;
    }
    double energy() const
    {
      return m_energy;
    }
    int detector_id() const
    {
      return m_detector_id;
    }

    inline void *operator new(size_t);
    inline void operator delete(void *hit);
  };

  using PhotocathodeHitsCollection = G4THitsCollection<PhotocathodeHit>;

  extern G4ThreadLocal G4Allocator<PhotocathodeHit> *PhotocathodeHitAllocator;

  inline void *PhotocathodeHit::operator new(size_t)
  {
    if (PhotocathodeHitAllocator == nullptr)
      PhotocathodeHitAllocator = new G4Allocator<PhotocathodeHit>;
    return static_cast<void *>(PhotocathodeHitAllocator->MallocSingle());
  }

  inline void PhotocathodeHit::operator delete(void *hit)
  {
    PhotocathodeHitAllocator->FreeSingle(static_cast<PhotocathodeHit *>(hit));
  }

} // namespace riptide

#endif
