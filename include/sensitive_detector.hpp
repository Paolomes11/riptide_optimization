#ifndef RIPTIDE_SENSITIVE_PHOTOCATHODE_HPP
#define RIPTIDE_SENSITIVE_PHOTOCATHODE_HPP

#include "photocathode_hit.hpp"

#include <G4VSensitiveDetector.hh>

namespace riptide
{

    class SensitivePhotocathode : public G4VSensitiveDetector
    {
        PhotocathodeHitsCollection *m_hits_collection{nullptr};
        int m_detector_id{0};
        int m_hc_id{-1};

    public:
        SensitivePhotocathode(G4String const &name, int detector_id);
        using G4VSensitiveDetector::G4VSensitiveDetector;
        
        void Initialize(G4HCofThisEvent *hce) override;
        bool ProcessHits(G4Step *step, G4TouchableHistory *history) override;
    };

} // namespace riptide

#endif