#include "psf_dof_action_initialization.hpp"

#include "primary_generator_action.hpp"
#include "psf_dof_event_action.hpp"
#include "psf_dof_stepping_action.hpp"

#include "optimization/detector_construction.hpp"

#include <G4RunManager.hh>

namespace riptide {

void PsfDofActionInitialization::BuildForMaster() const {
}

void PsfDofActionInitialization::Build() const {
  auto* primaryGen = new PrimaryGeneratorAction();
  if (m_useImportanceSampling) {
    auto* det = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    primaryGen->SetImportanceSampling(true, [det]() { return det->GetLens75Params(); });
  }
  SetUserAction(primaryGen);

  auto* eventAction = new PsfDofEventAction(false);
  SetUserAction(eventAction);

  SetUserAction(new PsfDofSteppingAction(eventAction));
}

} // namespace riptide
