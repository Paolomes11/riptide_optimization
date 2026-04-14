#ifndef RIPTIDE_DOF_STEPPING_ACTION_HPP
#define RIPTIDE_DOF_STEPPING_ACTION_HPP

#include <G4UserSteppingAction.hh>

namespace riptide {

class DofEventAction;

class DofSteppingAction : public G4UserSteppingAction {
 public:
  explicit DofSteppingAction(DofEventAction* event_action);
  ~DofSteppingAction() override = default;

  void SetLensAperturePlanes(double x1_mm, double r1_mm, double x2_mm, double r2_mm) {
    m_xLens1Aperture   = x1_mm;
    m_rLens1Aperture   = r1_mm;
    m_xLens2Aperture   = x2_mm;
    m_rLens2Aperture   = r2_mm;
    m_hasLensApertures = true;
  }

  void SetVirtualPlane(double x_mm) {
    m_xVirtual = x_mm;
  }

  void UserSteppingAction(const G4Step* step) override;

 private:
  double m_xVirtual             = 0.0;
  double m_xLens1Aperture       = 0.0;
  double m_rLens1Aperture       = 0.0;
  double m_xLens2Aperture       = 0.0;
  double m_rLens2Aperture       = 0.0;
  bool m_hasLensApertures       = false;
  bool m_crossedLens1           = false;
  DofEventAction* m_eventAction = nullptr;
};

} // namespace riptide

#endif // RIPTIDE_DOF_STEPPING_ACTION_HPP
