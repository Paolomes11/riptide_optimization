#ifndef RIPTIDE_PSF_DOF_ACTION_INITIALIZATION_HPP
#define RIPTIDE_PSF_DOF_ACTION_INITIALIZATION_HPP

#include <G4VUserActionInitialization.hh>

namespace riptide {

class PsfDofActionInitialization : public G4VUserActionInitialization {
 public:
  explicit PsfDofActionInitialization(bool use_importance_sampling)
      : G4VUserActionInitialization()
      , m_useImportanceSampling(use_importance_sampling) {
  }

  ~PsfDofActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build() const override;

 private:
  bool m_useImportanceSampling = false;
};

} // namespace riptide

#endif // RIPTIDE_PSF_DOF_ACTION_INITIALIZATION_HPP
