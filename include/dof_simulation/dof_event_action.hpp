#ifndef RIPTIDE_DOF_EVENT_ACTION_HPP
#define RIPTIDE_DOF_EVENT_ACTION_HPP

#include <G4UserEventAction.hh>

#include <vector>

namespace riptide {

class DofEventAction : public G4UserEventAction {
 public:
  DofEventAction()           = default;
  ~DofEventAction() override = default;

  void SetConfigId(int config_id) {
    m_configId = config_id;
  }

  int GetConfigId() const {
    return m_configId;
  }

  void AddRay(double y, double z, double dy, double dz, double w, double y_source);
  void ClearRays();

  double GetWeightedRayCount() const {
    return m_weightedRayCount;
  }

  std::vector<float>& YHits() {
    return m_yHits;
  }
  std::vector<float>& ZHits() {
    return m_zHits;
  }
  std::vector<float>& DyHits() {
    return m_dyHits;
  }
  std::vector<float>& DzHits() {
    return m_dzHits;
  }
  std::vector<float>& WeightHits() {
    return m_weightHits;
  }
  std::vector<float>& YSourceHits() {
    return m_ySourceHits;
  }

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  static DofEventAction* GetEventAction() {
    return s_currentEventAction;
  }

 private:
  int m_configId            = -1;
  double m_weightedRayCount = 0.0;

  std::vector<float> m_yHits;
  std::vector<float> m_zHits;
  std::vector<float> m_dyHits;
  std::vector<float> m_dzHits;
  std::vector<float> m_weightHits;
  std::vector<float> m_ySourceHits;

  static inline DofEventAction* s_currentEventAction = nullptr;
};

} // namespace riptide

#endif // RIPTIDE_DOF_EVENT_ACTION_HPP
