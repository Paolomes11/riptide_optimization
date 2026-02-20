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

    public:
        explicit DetectorConstruction(std::filesystem::path geometry_path);
        G4VPhysicalVolume *Construct() override;
        // void ConstructSDandField() override;
    };

} // namespace riptide

#endif // RIPTIDE_DETECTOR_CONSTRUCTION_HPP