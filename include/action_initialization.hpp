#ifndef RIPTIDE_ACTION_INITIALIZATION_HPP
#define RIPTIDE_ACTION_INITIALIZATION_HPP

#include <G4VUserActionInitialization.hh>

namespace riptide
{

    class ActionInitialization : public G4VUserActionInitialization
    {
    public:
        void BuildForMaster() const override;
        void Build() const override;
    };

} // namespace riptide

#endif // RIPTIDE_ACTION_INITIALIZATION_HPP