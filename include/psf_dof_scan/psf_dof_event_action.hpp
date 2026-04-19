#ifndef RIPTIDE_PSF_DOF_EVENT_ACTION_HPP
#define RIPTIDE_PSF_DOF_EVENT_ACTION_HPP

#include <G4UserEventAction.hh>

#include <vector>

namespace riptide {

struct PsfDofMoments {
  double sum_w     = 0.0;
  double sum_y     = 0.0;
  double sum_z     = 0.0;
  double sum_dy    = 0.0;
  double sum_dz    = 0.0;
  double sum_yy    = 0.0;
  double sum_zz    = 0.0;
  double sum_yz    = 0.0;
  double sum_dy_dy = 0.0;
  double sum_dz_dz = 0.0;
  double sum_y_dy  = 0.0;
  double sum_z_dz  = 0.0;
};

class PsfDofEventAction : public G4UserEventAction {
 public:
  explicit PsfDofEventAction(bool save_hits);
  ~PsfDofEventAction() override = default;

  void SetSaveHits(bool save_hits) {
    m_saveHits = save_hits;
  }

  void AddRay(double y, double z, double dy, double dz, double w);
  void ClearRays();

  const PsfDofMoments& Moments() const {
    return m_moments;
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

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

 private:
  bool m_saveHits = false;
  PsfDofMoments m_moments;

  std::vector<float> m_yHits;
  std::vector<float> m_zHits;
  std::vector<float> m_dyHits;
  std::vector<float> m_dzHits;
  std::vector<float> m_weightHits;
};

} // namespace riptide

#endif // RIPTIDE_PSF_DOF_EVENT_ACTION_HPP
