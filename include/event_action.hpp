#ifndef RIPTIDE_EVENT_ACTION_HPP
#define RIPTIDE_EVENT_ACTION_HPP

#include <G4UserEventAction.hh>

namespace riptide {

class EfficiencyCollector;

class EventAction : public G4UserEventAction {
  EfficiencyCollector* m_collector;

 public:
  EventAction(EfficiencyCollector* collector);

  void EndOfEventAction(const G4Event* event) override;
};

} // namespace riptide

#endif