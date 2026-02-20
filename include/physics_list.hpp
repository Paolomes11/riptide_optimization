#ifndef RIPTIDE_PHYSICS_LIST_HPP
#define RIPTIDE_PHYSICS_LIST_HPP

#include <G4VModularPhysicsList.hh>

namespace riptide
{

    class PhysicsList : public G4VModularPhysicsList
    {
    public:
        PhysicsList();
        ~PhysicsList() override = default;
    };

} // namespace riptide

#endif // RIPTIDE_PHYSICS_LIST_HPP