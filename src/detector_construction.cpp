#include "detector_construction.hpp"
#include "sensitive_detector.hpp"

#include <G4LogicalVolumeStore.hh>
#include <G4SDManager.hh>

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

    void DetectorConstruction::ConstructSDandField()
    {
        auto volume_store = G4LogicalVolumeStore::GetInstance();
        auto *sd_manager = G4SDManager::GetSDMpointer();

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

                // Estrai l'ID del detector dal nome (es. "photocathode_x" -> 0, "photocathode_y" -> 1)
                int detector_id = 0;
                if (it->value.find("_x") != std::string::npos)
                    detector_id = 0;
                else if (it->value.find("_y") != std::string::npos)
                    detector_id = 1;
                else if (it->value.find("_z") != std::string::npos)
                    detector_id = 2;

                // Crea il sensitive detector
                auto *sd = new SensitivePhotocathode{it->value, detector_id};
                // Registra il sensitive detector con il G4SDManager
                sd_manager->AddNewDetector(sd);
                // Associa il sensitive detector al logical volume
                lv->SetSensitiveDetector(sd);
            }
        }
    }

}