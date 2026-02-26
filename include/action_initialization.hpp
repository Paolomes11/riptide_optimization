#ifndef RIPTIDE_ACTION_INITIALIZATION_HPP
#define RIPTIDE_ACTION_INITIALIZATION_HPP

#include <G4VUserActionInitialization.hh>

namespace riptide {

class EfficiencyCollector;

class ActionInitialization : public G4VUserActionInitialization {
  EfficiencyCollector* m_collector;

 public:
  explicit ActionInitialization(EfficiencyCollector* collector = nullptr);

  void BuildForMaster() const override;
  void Build() const override;
};

} // namespace riptide

#endif // RIPTIDE_ACTION_INITIALIZATION_HPP