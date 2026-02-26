#include "efficiency_collector.hpp"

namespace riptide {

// Static member definition
EfficiencyCollector* EfficiencyCollector::s_instance = nullptr;

void EfficiencyCollector::recordEvent(int event_id) {
  m_events_with_hits.insert(event_id);
}

void EfficiencyCollector::recordEventWithHit(int event_id) {
  m_events_with_hits.insert(event_id);
}

double EfficiencyCollector::computeEfficiency(int n_photons_shot) const {
  if (n_photons_shot == 0)
    return 0.0;
  return static_cast<double>(m_events_with_hits.size()) / n_photons_shot;
}

void EfficiencyCollector::reset() {
  m_events_with_hits.clear();
}

// Global access methods
EfficiencyCollector* EfficiencyCollector::GetInstance() {
  return s_instance;
}

void EfficiencyCollector::SetInstance(EfficiencyCollector* collector) {
  s_instance = collector;
}

} // namespace riptide
