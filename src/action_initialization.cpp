#include "action_initialization.hpp"
#include "event_action.hpp"
#include "primary_generator_action.hpp"
#include "run_action.hpp"
#include <iostream>

namespace riptide {

ActionInitialization::ActionInitialization(const std::string& output_file)
    : m_output_file(output_file) {
  std::cout << "ActionInitialization: Constructor called, output file: " << output_file
            << std::endl;
}

void ActionInitialization::BuildForMaster() const {
  std::cout << "ActionInitialization: BuildForMaster called" << std::endl;
}

void ActionInitialization::Build() const {
  // Azioni standard
  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new RunAction());

  // EventAction con file ROOT
  SetUserAction(new EventAction(m_output_file));
}

} // namespace riptide
