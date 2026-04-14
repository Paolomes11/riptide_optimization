#ifndef RIPTIDE_DOF_ACTION_INITIALIZATION_HPP
#define RIPTIDE_DOF_ACTION_INITIALIZATION_HPP

#include <G4VUserActionInitialization.hh>

namespace riptide {

class DofActionInitialization : public G4VUserActionInitialization {
 public:
  explicit DofActionInitialization(bool use_importance_sampling)
      : G4VUserActionInitialization()
      , m_useImportanceSampling(use_importance_sampling) {
  }

  ~DofActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build() const override;

 private:
  bool m_useImportanceSampling = false;
};

} // namespace riptide

#endif // RIPTIDE_DOF_ACTION_INITIALIZATION_HPP
