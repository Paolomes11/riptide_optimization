#include "event_action.hpp"
#include "efficiency_collector.hpp"
#include "G4Event.hh"
#include <iostream>

namespace riptide {

EventAction::EventAction(riptide::EfficiencyCollector* collector)
    : m_collector(collector) {
}

void EventAction::EndOfEventAction(const G4Event* event) {

}

} // namespace riptide
