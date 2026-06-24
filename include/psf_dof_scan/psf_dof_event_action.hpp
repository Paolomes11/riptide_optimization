#ifndef RIPTIDE_PSF_DOF_EVENT_ACTION_HPP
#define RIPTIDE_PSF_DOF_EVENT_ACTION_HPP

#include "spot_grid.hpp"

#include <G4UserEventAction.hh>

#include <cmath>
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

  // Legacy (un run = uno spot)
  void AddRay(double y, double z, double dy, double dz, double w);
  void ClearRays();

  const PsfDofMoments& Moments() const { return m_moments; }

  std::vector<float>& YHits()      { return m_yHits; }
  std::vector<float>& ZHits()      { return m_zHits; }
  std::vector<float>& DyHits()     { return m_dyHits; }
  std::vector<float>& DzHits()     { return m_dzHits; }
  std::vector<float>& WeightHits() { return m_weightHits; }

  // Spot mode (un BeamOn per coppia, tutti gli spot in parallelo)
  bool IsSpotMode() const { return m_spotMode; }

  int CalcSpotId(double x_src_mm, double y_src_mm) const {
    return calc_spot_id(x_src_mm, y_src_mm, m_xSrcMin, m_dxSrc, m_ySrcMin, m_dySrc, m_nySrc);
  }

  void InitSpotMode(int n_spots, double x_src_min, double dx_src,
                    double y_src_min, double dy_src, int ny_src);
  void ResetSpotAccumulators();
  void AddRay(int spot_id, double y, double z, double dy, double dz, double w);
  void AddKilledLens1(int spot_id);
  void AddKilledLens2(int spot_id);
  void AddKilledBack(int spot_id);

  const std::vector<PsfDofMoments>& GetSpotMoments()     const { return m_spotMoments; }
  const std::vector<int>& GetSpotKilledLens1()           const { return m_spotKilledLens1; }
  const std::vector<int>& GetSpotKilledLens2()           const { return m_spotKilledLens2; }
  const std::vector<int>& GetSpotKilledBack()            const { return m_spotKilledBack; }

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

  // Spot mode accumulators
  bool m_spotMode = false;
  std::vector<PsfDofMoments> m_spotMoments;
  std::vector<int> m_spotKilledLens1;
  std::vector<int> m_spotKilledLens2;
  std::vector<int> m_spotKilledBack;
  double m_xSrcMin = 0.0;
  double m_dxSrc   = 1.0;
  double m_ySrcMin = 0.0;
  double m_dySrc   = 1.0;
  int    m_nySrc   = 1;
};

} // namespace riptide

#endif // RIPTIDE_PSF_DOF_EVENT_ACTION_HPP
