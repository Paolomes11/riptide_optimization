#include "physics_list.hpp"

#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4HadronPhysicsINCLXX.hh>
#include <G4OpticalParameters.hh>
#include <G4OpticalPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4SystemOfUnits.hh>

namespace riptide
{

    PhysicsList::PhysicsList()
    {
        G4VModularPhysicsList::defaultCutValue = 1.0 * mm;

        RegisterPhysics(new G4EmStandardPhysics{0});
        RegisterPhysics(new G4HadronElasticPhysics{0});
        RegisterPhysics(
            new G4HadronPhysicsINCLXX("hadronPhysicsINCLXX"));
        RegisterPhysics(new G4DecayPhysics{0});
        RegisterPhysics(new G4StoppingPhysics{0});
        RegisterPhysics(new G4OpticalPhysics{0});
    }

} // namespace riptide
