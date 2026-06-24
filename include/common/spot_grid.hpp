#ifndef RIPTIDE_SPOT_GRID_HPP
#define RIPTIDE_SPOT_GRID_HPP

#include "primary_generator_action.hpp"

#include <cmath>
#include <vector>

namespace riptide {

inline int calc_spot_id(double x_mm, double y_mm, double x_min, double dx, double y_min,
                        double dy, int n_y) {
  int ix = static_cast<int>(std::lround((x_mm - x_min) / dx));
  int iy = static_cast<int>(std::lround((y_mm - y_min) / dy));
  return ix * n_y + iy;
}

inline std::vector<SpotConfig> build_spot_grid(double x_min, double x_max, double dx,
                                               double y_min, double y_max, double dy) {
  std::vector<SpotConfig> spots;
  for (double x = x_min; x <= x_max + 1e-9; x += dx) {
    for (double y = y_min; y <= y_max + 1e-9; y += dy) {
      SpotConfig sc;
      sc.pos    = G4ThreeVector(x, y, 0.0);
      sc.use_is = false;
      spots.push_back(sc);
    }
  }
  return spots;
}

} // namespace riptide

#endif // RIPTIDE_SPOT_GRID_HPP
