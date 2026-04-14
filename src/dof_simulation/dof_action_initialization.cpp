#include "dof_action_initialization.hpp"

#include "dof_event_action.hpp"
#include "dof_stepping_action.hpp"
#include "primary_generator_action.hpp"

#include "optimization/detector_construction.hpp"

#include <G4RunManager.hh>

namespace riptide {

void DofActionInitialization::BuildForMaster() const {
}

void DofActionInitialization::Build() const {
  auto* primaryGen = new PrimaryGeneratorAction();
  if (m_useImportanceSampling) {
    auto* det = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    primaryGen->SetImportanceSampling(true, [det]() { return det->GetLens75Params(); });
  }
  SetUserAction(primaryGen);

  auto* eventAction = new DofEventAction();
  SetUserAction(eventAction);

  SetUserAction(new DofSteppingAction(eventAction));
}

} // namespace riptide
