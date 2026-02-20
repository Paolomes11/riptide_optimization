#include "detector_construction.hpp"
// #include "sensitive_photocathode.hpp"

#include <G4LogicalVolumeStore.hh>

namespace riptide
{

    DetectorConstruction::DetectorConstruction(std::filesystem::path geometry_path)
        : m_geometry_path{std::move(geometry_path)}
    {
        if (!std::filesystem::exists(m_geometry_path) || !std::filesystem::is_regular_file(m_geometry_path))
        {
            throw std::runtime_error("Geometry path is invalid");
        }
    }

    G4VPhysicalVolume *DetectorConstruction::Construct()
    {
        m_parser.Read(m_geometry_path.string());
        return m_parser.GetWorldVolume();
    }

    /*
    void DetectorConstruction::ConstructSDandField()
    {
        auto volume_store = G4LogicalVolumeStore::GetInstance();

        for (auto lv : *volume_store)
        {
            auto const &aux_list = m_parser.GetVolumeAuxiliaryInformation(lv);
            auto const it = std::find_if(aux_list.begin(), aux_list.end(),
                                         [](const auto &aux)
                                         { return aux.type == "sd_name"; });
            if (it != aux_list.end())
            {
                if (it->value.empty())
                {
                    throw std::runtime_error("Empty sd_name value for volume");
                }

                lv->SetSensitiveDetector(new SensitivePhotocathode{it->value});
            }
        }
    }*/

}