#ifndef RIPTIDE_SENSITIVE_PHOTOCATHODE_HPP
#define RIPTIDE_SENSITIVE_PHOTOCATHODE_HPP

#include <G4VSensitiveDetector.hh>

namespace riptide
{

    class SensitivePhotocathode : public G4VSensitiveDetector
    {
    public:
        using G4VSensitiveDetector::G4VSensitiveDetector;
        bool ProcessHits(G4Step *step, G4TouchableHistory *history) override;
    };

} // namespace riptide

#endif