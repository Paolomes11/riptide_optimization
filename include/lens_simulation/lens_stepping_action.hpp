#ifndef RIPTIDE_LENS_STEPPING_ACTION_HPP
#define RIPTIDE_LENS_STEPPING_ACTION_HPP

#include <G4UserSteppingAction.hh>

namespace riptide {

class EventAction;

class LensSteppingAction : public G4UserSteppingAction {
 public:
  explicit LensSteppingAction(EventAction* event_action);
  ~LensSteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

  void SetVirtualPlane(double x_mm) { m_xVirtual = x_mm; }

  void SetLensAperturePlanes(double x1_mm, double r1_mm, double x2_mm, double r2_mm) {
    m_xLens1Aperture   = x1_mm;
    m_rLens1Aperture   = r1_mm;
    m_xLens2Aperture   = x2_mm;
    m_rLens2Aperture   = r2_mm;
    m_hasLensApertures = true;
  }

 private:
  EventAction* m_eventAction   = nullptr;
  double m_xVirtual            = 1e9;
  bool   m_hasLensApertures    = false;
  double m_xLens1Aperture      = 0.0;
  double m_rLens1Aperture      = 0.0;
  double m_xLens2Aperture      = 0.0;
  double m_rLens2Aperture      = 0.0;
  bool   m_crossedLens1        = false;
};

} // namespace riptide

#endif // RIPTIDE_LENS_STEPPING_ACTION_HPP
