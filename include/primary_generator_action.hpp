#ifndef RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP
#define RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP

#include <G4VUserPrimaryGeneratorAction.hh>

// Forward declarations for compilation speedup
class G4Event;
class G4GeneralParticleSource;

namespace riptide
{

    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
        // m_prefix is a member variable
        G4GeneralParticleSource *m_gps{nullptr};

    public:
        PrimaryGeneratorAction();
        ~PrimaryGeneratorAction() override;

        void GeneratePrimaries(G4Event *event) override;
    };

} // namespace riptide

#endif // RIPTIDE_PRIMARY_GENERATOR_ACTION_HPP