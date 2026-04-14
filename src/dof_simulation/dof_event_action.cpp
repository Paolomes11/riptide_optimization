#include "dof_event_action.hpp"

#include <G4Event.hh>

namespace riptide {

void DofEventAction::AddRay(double y, double z, double dy, double dz, double w, double y_source) {
  m_yHits.push_back(static_cast<float>(y));
  m_zHits.push_back(static_cast<float>(z));
  m_dyHits.push_back(static_cast<float>(dy));
  m_dzHits.push_back(static_cast<float>(dz));
  m_weightHits.push_back(static_cast<float>(w));
  m_ySourceHits.push_back(static_cast<float>(y_source));
  m_weightedRayCount += w;
}

void DofEventAction::ClearRays() {
  m_weightedRayCount = 0.0;
  m_yHits.clear();
  m_zHits.clear();
  m_dyHits.clear();
  m_dzHits.clear();
  m_weightHits.clear();
  m_ySourceHits.clear();
}

void DofEventAction::BeginOfEventAction(const G4Event* /*event*/) {
  s_currentEventAction = this;
}

void DofEventAction::EndOfEventAction(const G4Event* /*event*/) {
}

} // namespace riptide
