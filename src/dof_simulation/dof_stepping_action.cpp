#include "dof_stepping_action.hpp"

#include "dof_event_action.hpp"

#include <G4OpticalPhoton.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>

#include <cmath>

namespace riptide {

DofSteppingAction::DofSteppingAction(DofEventAction* event_action)
    : G4UserSteppingAction()
    , m_eventAction(event_action) {
}

void DofSteppingAction::UserSteppingAction(const G4Step* step) {
  auto* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
    return;
  }

  auto pos      = track->GetPosition();
  auto momentum = track->GetMomentum();

  if (pos.x() < -1.0 * CLHEP::mm && momentum.x() < 0.0) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  if (std::abs(pos.y()) > 150.0 * CLHEP::mm || std::abs(pos.z()) > 150.0 * CLHEP::mm) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  double pre_x  = step->GetPreStepPoint()->GetPosition().x() / CLHEP::mm;
  double post_x = step->GetPostStepPoint()->GetPosition().x() / CLHEP::mm;

  if (pre_x >= m_xVirtual || post_x < m_xVirtual) {
    return;
  }

  auto pp = step->GetPreStepPoint()->GetPosition();
  auto qp = step->GetPostStepPoint()->GetPosition();

  double denom = (post_x - pre_x);
  if (std::abs(denom) < 1e-12) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  double t = (m_xVirtual - pre_x) / denom;

  double y = (pp.y() + t * (qp.y() - pp.y())) / CLHEP::mm;
  double z = (pp.z() + t * (qp.z() - pp.z())) / CLHEP::mm;

  auto mom_dir = track->GetMomentumDirection();
  if (std::abs(mom_dir.x()) < 1e-6) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  double dy       = mom_dir.y() / mom_dir.x();
  double dz       = mom_dir.z() / mom_dir.x();
  double w        = track->GetWeight();
  double y_source = track->GetVertexPosition().y() / CLHEP::mm;

  if (m_eventAction) {
    m_eventAction->AddRay(y, z, dy, dz, w, y_source);
  }

  track->SetTrackStatus(fStopAndKill);
}

} // namespace riptide
