#ifndef RIPTIDE_DETECTOR_CONSTRUCTION_HPP
#define RIPTIDE_DETECTOR_CONSTRUCTION_HPP

#include <G4VUserDetectorConstruction.hh>
#include <G4GDMLParser.hh>

#include <filesystem>

namespace riptide
{

    class DetectorConstruction : public G4VUserDetectorConstruction
    {
        G4GDMLParser m_parser{};
        std::filesystem::path m_geometry_path;

        // Posizioni delle lenti (per ottimizzazione)
        double m_lens75_x;
        double m_lens60_x;

    public:
        DetectorConstruction(std::filesystem::path geometry_path, double lens75_x = 83.9, double lens60_x = 153.4);
        G4VPhysicalVolume *Construct() override;
        void ConstructSDandField() override;
    };

} // namespace riptide

#endif // RIPTIDE_DETECTOR_CONSTRUCTION_HPP