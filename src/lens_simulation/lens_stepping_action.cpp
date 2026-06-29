#include "lens_stepping_action.hpp"
#include "event_action.hpp"

#include <G4OpticalPhoton.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>

#include <cmath>

namespace riptide {

LensSteppingAction::LensSteppingAction(EventAction* event_action)
    : G4UserSteppingAction()
    , m_eventAction(event_action) {}

void LensSteppingAction::UserSteppingAction(const G4Step* step) {
  auto* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
    return;
  }

  // Reset flag per questa traccia (lo step 1 = primo step del fotone)
  if (track->GetCurrentStepNumber() <= 1) {
    m_crossedLens1 = false;
  }

  // Calcola spot_id dalla posizione del vertice primario
  int spot_id = -1;
  if (m_eventAction && m_eventAction->IsSpotMode()) {
    G4ThreeVector vtx = track->GetVertexPosition();
    spot_id = m_eventAction->CalcSpotId(vtx.x() / CLHEP::mm, vtx.y() / CLHEP::mm);
  }

  // ── Kill verso la sorgente (backward) ─────────────────────────────────
  auto pos      = track->GetPosition();
  auto momentum = track->GetMomentum();
  if (pos.x() < -1.0 * CLHEP::mm && momentum.x() < 0.0) {
    if (m_eventAction) m_eventAction->AddKilledVBack(spot_id);
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  double pre_x  = step->GetPreStepPoint()->GetPosition().x() / CLHEP::mm;
  double post_x = step->GetPostStepPoint()->GetPosition().x() / CLHEP::mm;

  // ── Kill per diaframma L1 ──────────────────────────────────────────────
  if (m_hasLensApertures && pre_x < m_xLens1Aperture && post_x >= m_xLens1Aperture) {
    auto pp      = step->GetPreStepPoint()->GetPosition();
    auto qp      = step->GetPostStepPoint()->GetPosition();
    double denom = post_x - pre_x;
    if (std::abs(denom) < 1e-12) {
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    double t = (m_xLens1Aperture - pre_x) / denom;
    double y = (pp.y() + t * (qp.y() - pp.y())) / CLHEP::mm;
    double z = (pp.z() + t * (qp.z() - pp.z())) / CLHEP::mm;
    if (std::hypot(y, z) > m_rLens1Aperture) {
      if (m_eventAction) m_eventAction->AddKilledVLens1(spot_id);
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    m_crossedLens1 = true;
  }

  // ── Kill per diaframma L2 ──────────────────────────────────────────────
  if (m_hasLensApertures && m_crossedLens1 && pre_x < m_xLens2Aperture
      && post_x >= m_xLens2Aperture) {
    auto pp      = step->GetPreStepPoint()->GetPosition();
    auto qp      = step->GetPostStepPoint()->GetPosition();
    double denom = post_x - pre_x;
    if (std::abs(denom) < 1e-12) {
      track->SetTrackStatus(fStopAndKill);
      return;
    }
    double t = (m_xLens2Aperture - pre_x) / denom;
    double y = (pp.y() + t * (qp.y() - pp.y())) / CLHEP::mm;
    double z = (pp.z() + t * (qp.z() - pp.z())) / CLHEP::mm;
    if (std::hypot(y, z) > m_rLens2Aperture) {
      if (m_eventAction) m_eventAction->AddKilledVLens2(spot_id);
      track->SetTrackStatus(fStopAndKill);
      return;
    }
  }

  // ── Registrazione al piano virtuale ───────────────────────────────────
  // Procede solo quando il fotone attraversa x_virtual da sinistra (una volta sola)
  if (pre_x >= m_xVirtual || post_x < m_xVirtual) return;

  auto pp      = step->GetPreStepPoint()->GetPosition();
  auto qp      = step->GetPostStepPoint()->GetPosition();
  double denom = post_x - pre_x;
  if (std::abs(denom) < 1e-12) return;

  double t = (m_xVirtual - pre_x) / denom;
  double y = (pp.y() + t * (qp.y() - pp.y())) / CLHEP::mm;
  double z = (pp.z() + t * (qp.z() - pp.z())) / CLHEP::mm;

  auto mom_dir = track->GetMomentumDirection();
  if (std::abs(mom_dir.x()) < 1e-6) return;
  double dy = mom_dir.y() / mom_dir.x();
  double dz = mom_dir.z() / mom_dir.x();
  double w  = track->GetWeight();

  if (m_eventAction)
    m_eventAction->AddVirtualHit(spot_id, y, z, dy, dz, w);

  // NON uccidiamo il fotone — continua verso il fotocatodo
}

} // namespace riptide
