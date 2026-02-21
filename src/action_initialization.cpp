#include "action_initialization.hpp"
#include "primary_generator_action.hpp"
#include "run_action.hpp"
#include "event_action.hpp"

namespace riptide
{

    void ActionInitialization::BuildForMaster() const
    {
        
    }

    void ActionInitialization::Build() const
    {
        // Initialize user actions for worker threads here
        // For example:
        SetUserAction(new PrimaryGeneratorAction());
        SetUserAction(new RunAction());
        SetUserAction(new EventAction());
        // SetUserAction(new SteppingAction());
    }

} // namespace riptide