#include "psf_dof_stepping_action.hpp"

#include "psf_dof_event_action.hpp"

#include <G4OpticalPhoton.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>

#include <cmath>

namespace riptide {

PsfDofSteppingAction::PsfDofSteppingAction(PsfDofEventAction* event_action)
    : G4UserSteppingAction()
    , m_eventAction(event_action) {
}

void PsfDofSteppingAction::UserSteppingAction(const G4Step* step) {
  auto* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
    return;
  }

  if (track->GetCurrentStepNumber() <= 1) {
    m_crossedLens1 = false;
  }

  auto pos      = track->GetPosition();
  auto momentum = track->GetMomentum();

  // Calcola spot_id una sola volta per questa traccia (O(1) via vertex position)
  int spot_id = -1;
  if (m_eventAction && m_eventAction->IsSpotMode()) {
    G4ThreeVector vtx = track->GetVertexPosition();
    spot_id = m_eventAction->CalcSpotId(vtx.x() / CLHEP::mm, vtx.y() / CLHEP::mm);
  }

  if (pos.x() < -1.0 * CLHEP::mm && momentum.x() < 0.0) {
    if (spot_id >= 0) {
      m_eventAction->AddKilledBack(spot_id);
    } else {
      ++m_n_killed_back;
    }
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  double pre_x  = step->GetPreStepPoint()->GetPosition().x() / CLHEP::mm;
  double post_x = step->GetPostStepPoint()->GetPosition().x() / CLHEP::mm;

  if (m_hasLensApertures && pre_x < m_xLens1Aperture && post_x >= m_xLens1Aperture) {
    auto pp      = step->GetPreStepPoint()->GetPosition();
    auto qp      = step->GetPostStepPoint()->GetPosition();
    double denom = (post_x - pre_x);
    if (std::abs(denom) < 1e-12) {
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    double t = (m_xLens1Aperture - pre_x) / denom;
    double y = (pp.y() + t * (qp.y() - pp.y())) / CLHEP::mm;
    double z = (pp.z() + t * (qp.z() - pp.z())) / CLHEP::mm;
    double r = std::hypot(y, z);
    if (r > m_rLens1Aperture) {
      if (spot_id >= 0) {
        m_eventAction->AddKilledLens1(spot_id);
      } else {
        ++m_n_killed_lens1;
      }
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    m_crossedLens1 = true;
  }

  if (m_hasLensApertures && m_crossedLens1 && pre_x < m_xLens2Aperture
      && post_x >= m_xLens2Aperture) {
    auto pp      = step->GetPreStepPoint()->GetPosition();
    auto qp      = step->GetPostStepPoint()->GetPosition();
    double denom = (post_x - pre_x);
    if (std::abs(denom) < 1e-12) {
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    double t = (m_xLens2Aperture - pre_x) / denom;
    double y = (pp.y() + t * (qp.y() - pp.y())) / CLHEP::mm;
    double z = (pp.z() + t * (qp.z() - pp.z())) / CLHEP::mm;
    double r = std::hypot(y, z);
    if (r > m_rLens2Aperture) {
      if (spot_id >= 0) {
        m_eventAction->AddKilledLens2(spot_id);
      } else {
        ++m_n_killed_lens2;
      }
      track->SetTrackStatus(fStopAndKill);
      return;
    }
  }

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

  double dy = mom_dir.y() / mom_dir.x();
  double dz = mom_dir.z() / mom_dir.x();
  double w  = track->GetWeight();

  if (m_eventAction) {
    if (spot_id >= 0) {
      m_eventAction->AddRay(spot_id, y, z, dy, dz, w);
    } else {
      m_eventAction->AddRay(y, z, dy, dz, w);
    }
  }

  track->SetTrackStatus(fStopAndKill);
}

} // namespace riptide
