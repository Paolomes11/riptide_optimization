#include "detector_construction.hpp"
#include "sensitive_detector.hpp"

#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>

namespace riptide
{

    DetectorConstruction::DetectorConstruction(std::filesystem::path geometry_path, double lens75_x, double lens60_x)
        : m_geometry_path{std::move(geometry_path)}, m_lens75_x{lens75_x}, m_lens60_x{lens60_x}
    {
        if (!std::filesystem::exists(m_geometry_path) || !std::filesystem::is_regular_file(m_geometry_path))
        {
            throw std::runtime_error("Geometry path is invalid");
        }
    }

    G4VPhysicalVolume *DetectorConstruction::Construct()
    {
        m_parser.Read(m_geometry_path.string());
        auto world = m_parser.GetWorldVolume();

        // Sposta le lenti per ottimizzazione
        auto pv_store = G4PhysicalVolumeStore::GetInstance();

        for (auto pv : *pv_store)
        {
            auto name = pv->GetName();
            if (name == "lens75_x_phys")
            {
                pv->SetTranslation(G4ThreeVector(m_lens75_x, 0., 0.));
            }
            else if (name == "lens60_x_phys")
            {
                pv->SetTranslation(G4ThreeVector(m_lens60_x, 0., 0.));
            }
        }

        return world;
    }

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
    }

}