#include "primary_generator_action.hpp"

#include <G4GeneralParticleSource.hh>

namespace riptide
{

    PrimaryGeneratorAction::PrimaryGeneratorAction()
        : m_gps{new G4GeneralParticleSource()}
    {
    }

    PrimaryGeneratorAction::~PrimaryGeneratorAction()
    {
        delete m_gps;
    }

    void PrimaryGeneratorAction::GeneratePrimaries(G4Event *event)
    {
        // Genera il primo vertice dell'interazione particella
        m_gps->GeneratePrimaryVertex(event);
    }

} // namespace riptide