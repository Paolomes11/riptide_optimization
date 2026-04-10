#ifndef RIPTIDE_IMPORTANCE_SAMPLING_HPP
#define RIPTIDE_IMPORTANCE_SAMPLING_HPP

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <algorithm>
#include <vector>

namespace riptide {

/**
 * @brief Helper class to calculate the optimal cone for geometric importance sampling.
 */
class ImportanceSamplingHelper {
 public:
  struct LensParams {
    double x;         // Position of the lens center (X axis)
    double diameter;  // Physical diameter
    double R1;        // Radius of curvature of the first surface (facing -X)
    double R2;        // Radius of curvature of the second surface (facing +X)
    double tc;        // Center thickness
    double te;        // Edge thickness
    bool is_biconvex; // True if biconvex, false if plano-convex
  };

  /**
   * @brief Calculates the cone axis and max opening angle to cover the lens from a set of source
   * points.
   *
   * @param sourcePoints List of source points (e.g. corners of a rectangle).
   * @param lens Lens parameters.
   * @param outAxis Output: axis of the cone (from average source position to lens center).
   * @param outMaxTheta Output: maximum half-opening angle of the cone.
   */
  static void CalculateGlobalCone(const std::vector<G4ThreeVector>& sourcePoints,
                                  const LensParams& lens, G4ThreeVector& outAxis,
                                  double& outMaxTheta);

  /**
   * @brief Calculates the max opening angle to cover the lens from a SINGLE source point,
   * relative to a PROVIDED axis.
   */
  static void CalculateCone(const G4ThreeVector& sourcePos, const LensParams& lens,
                            const G4ThreeVector& axis, double& outMaxTheta);

 private:
  static void checkSphereTangency(const G4ThreeVector& S, const G4ThreeVector& axis, double xc,
                                  double R, double r_lens, double& maxTheta);
  static void checkCylinderTangency(const G4ThreeVector& S, const G4ThreeVector& axis, double r,
                                    double x_min, double x_max, double& maxTheta);
};

} // namespace riptide

#endif // RIPTIDE_IMPORTANCE_SAMPLING_HPP
