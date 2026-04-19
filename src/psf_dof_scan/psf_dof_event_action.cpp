#include "psf_dof_event_action.hpp"

#include <G4Event.hh>

namespace riptide {

PsfDofEventAction::PsfDofEventAction(bool save_hits)
    : G4UserEventAction()
    , m_saveHits(save_hits) {
}

void PsfDofEventAction::AddRay(double y, double z, double dy, double dz, double w) {
  m_moments.sum_w += w;
  m_moments.sum_y += w * y;
  m_moments.sum_z += w * z;
  m_moments.sum_dy += w * dy;
  m_moments.sum_dz += w * dz;

  m_moments.sum_yy += w * y * y;
  m_moments.sum_zz += w * z * z;
  m_moments.sum_yz += w * y * z;

  m_moments.sum_dy_dy += w * dy * dy;
  m_moments.sum_dz_dz += w * dz * dz;

  m_moments.sum_y_dy += w * y * dy;
  m_moments.sum_z_dz += w * z * dz;

  if (m_saveHits) {
    m_yHits.push_back(static_cast<float>(y));
    m_zHits.push_back(static_cast<float>(z));
    m_dyHits.push_back(static_cast<float>(dy));
    m_dzHits.push_back(static_cast<float>(dz));
    m_weightHits.push_back(static_cast<float>(w));
  }
}

void PsfDofEventAction::ClearRays() {
  m_moments = {};
  m_yHits.clear();
  m_zHits.clear();
  m_dyHits.clear();
  m_dzHits.clear();
  m_weightHits.clear();
}

void PsfDofEventAction::BeginOfEventAction(const G4Event* /*event*/) {
}

void PsfDofEventAction::EndOfEventAction(const G4Event* /*event*/) {
}

} // namespace riptide
