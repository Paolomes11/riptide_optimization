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

void PsfDofEventAction::InitSpotMode(int n_spots, double x_src_min, double dx_src,
                                     double y_src_min, double dy_src, int ny_src) {
  m_spotMode = true;
  m_xSrcMin  = x_src_min;
  m_dxSrc    = dx_src;
  m_ySrcMin  = y_src_min;
  m_dySrc    = dy_src;
  m_nySrc    = ny_src;
  m_spotMoments.assign(n_spots, PsfDofMoments{});
  m_spotKilledLens1.assign(n_spots, 0);
  m_spotKilledLens2.assign(n_spots, 0);
  m_spotKilledBack.assign(n_spots, 0);
}

void PsfDofEventAction::ResetSpotAccumulators() {
  m_spotMode = false;
  m_spotMoments.clear();
  m_spotKilledLens1.clear();
  m_spotKilledLens2.clear();
  m_spotKilledBack.clear();
}

void PsfDofEventAction::AddRay(int spot_id, double y, double z, double dy, double dz, double w) {
  if (spot_id < 0 || spot_id >= static_cast<int>(m_spotMoments.size())) return;
  auto& m      = m_spotMoments[spot_id];
  m.sum_w     += w;
  m.sum_y     += w * y;
  m.sum_z     += w * z;
  m.sum_dy    += w * dy;
  m.sum_dz    += w * dz;
  m.sum_yy    += w * y  * y;
  m.sum_zz    += w * z  * z;
  m.sum_yz    += w * y  * z;
  m.sum_dy_dy += w * dy * dy;
  m.sum_dz_dz += w * dz * dz;
  m.sum_y_dy  += w * y  * dy;
  m.sum_z_dz  += w * z  * dz;
}

void PsfDofEventAction::AddKilledLens1(int spot_id) {
  if (spot_id >= 0 && spot_id < static_cast<int>(m_spotKilledLens1.size()))
    ++m_spotKilledLens1[spot_id];
}

void PsfDofEventAction::AddKilledLens2(int spot_id) {
  if (spot_id >= 0 && spot_id < static_cast<int>(m_spotKilledLens2.size()))
    ++m_spotKilledLens2[spot_id];
}

void PsfDofEventAction::AddKilledBack(int spot_id) {
  if (spot_id >= 0 && spot_id < static_cast<int>(m_spotKilledBack.size()))
    ++m_spotKilledBack[spot_id];
}

void PsfDofEventAction::BeginOfEventAction(const G4Event* /*event*/) {
}

void PsfDofEventAction::EndOfEventAction(const G4Event* /*event*/) {
}

} // namespace riptide
