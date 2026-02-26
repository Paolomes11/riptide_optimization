#include "action_initialization.hpp"
#include "efficiency_collector.hpp"
#include "event_action.hpp"
#include "primary_generator_action.hpp"
#include "run_action.hpp"
#include <iostream>

namespace riptide {

ActionInitialization::ActionInitialization(EfficiencyCollector* collector)
    : m_collector(collector) {
  std::cout << "ActionInitialization: Constructor called with collector: " << collector
            << std::endl;
}

void ActionInitialization::BuildForMaster() const {
  std::cout << "ActionInitialization: BuildForMaster called" << std::endl;
}

void ActionInitialization::Build() const {
  // Initialize user actions for worker threads here
  // For example:
  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new RunAction());
  // SetUserAction(new SteppingAction());

  // EventAction con collector, solo se presente
  if (m_collector) {
    SetUserAction(new EventAction(m_collector));
  }
}

} // namespace riptide
