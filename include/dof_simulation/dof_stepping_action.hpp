#ifndef RIPTIDE_DOF_STEPPING_ACTION_HPP
#define RIPTIDE_DOF_STEPPING_ACTION_HPP

#include <G4UserSteppingAction.hh>

namespace riptide {

class DofEventAction;

class DofSteppingAction : public G4UserSteppingAction {
 public:
  explicit DofSteppingAction(DofEventAction* event_action);
  ~DofSteppingAction() override = default;

  void SetVirtualPlane(double x_mm) {
    m_xVirtual = x_mm;
  }

  void UserSteppingAction(const G4Step* step) override;

 private:
  double m_xVirtual             = 0.0;
  DofEventAction* m_eventAction = nullptr;
};

} // namespace riptide

#endif // RIPTIDE_DOF_STEPPING_ACTION_HPP
