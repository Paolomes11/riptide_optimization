#ifndef RIPTIDE_PSF_DOF_STEPPING_ACTION_HPP
#define RIPTIDE_PSF_DOF_STEPPING_ACTION_HPP

#include <G4UserSteppingAction.hh>

namespace riptide {

class PsfDofEventAction;

class PsfDofSteppingAction : public G4UserSteppingAction {
 public:
  explicit PsfDofSteppingAction(PsfDofEventAction* event_action);
  ~PsfDofSteppingAction() override = default;

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
  double m_xVirtual                = 0.0;
  double m_xLens1Aperture          = 0.0;
  double m_rLens1Aperture          = 0.0;
  double m_xLens2Aperture          = 0.0;
  double m_rLens2Aperture          = 0.0;
  bool m_hasLensApertures          = false;
  bool m_crossedLens1              = false;
  PsfDofEventAction* m_eventAction = nullptr;
};

} // namespace riptide

#endif // RIPTIDE_PSF_DOF_STEPPING_ACTION_HPP
