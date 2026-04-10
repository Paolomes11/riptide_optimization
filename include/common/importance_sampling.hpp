#ifndef RIPTIDE_IMPORTANCE_SAMPLING_HPP
#define RIPTIDE_IMPORTANCE_SAMPLING_HPP

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <algorithm>
#include <cmath>
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
                                  double& outMaxTheta) {
    if (sourcePoints.empty())
      return;

    // Calculate average source position for the cone axis
    G4ThreeVector avgSource(0, 0, 0);
    for (const auto& p : sourcePoints)
      avgSource += p;
    avgSource /= static_cast<double>(sourcePoints.size());

    G4ThreeVector lensCenter(lens.x, 0, 0);
    outAxis = (lensCenter - avgSource).unit();

    double maxTheta = 0.0;
    for (const auto& s : sourcePoints) {
      double currentMaxTheta = 0.0;
      CalculateCone(s, lens, outAxis, currentMaxTheta);
      if (currentMaxTheta > maxTheta)
        maxTheta = currentMaxTheta;
    }
    outMaxTheta = maxTheta;
  }

  /**
   * @brief Calculates the max opening angle to cover the lens from a SINGLE source point,
   * relative to a PROVIDED axis.
   */
  static void CalculateCone(const G4ThreeVector& sourcePos, const LensParams& lens,
                            const G4ThreeVector& axis, double& outMaxTheta) {
    // We find the max angle by sampling critical points on the lens
    double maxTheta = 0.0;

    auto updateMaxTheta = [&](const G4ThreeVector& p) {
      G4ThreeVector dir = (p - sourcePos).unit();
      double cosTheta   = dir.dot(axis);
      cosTheta          = std::max(-1.0, std::min(1.0, cosTheta));
      double theta      = std::acos(cosTheta);
      if (theta > maxTheta)
        maxTheta = theta;
    };

    // 1. Sample the circular edges (where cylinder meets caps)
    const int nPhi = 36; // Sufficient for global cone
    double r       = lens.diameter / 2.0;

    // Use edge thickness for sampling the circular boundary
    double x_front_edge = lens.x - lens.te / 2.0;
    double x_back_edge  = lens.x + lens.te / 2.0;

    for (int i = 0; i < nPhi; ++i) {
      double phi = i * 2.0 * M_PI / nPhi;
      double y   = r * std::cos(phi);
      double z   = r * std::sin(phi);
      updateMaxTheta(G4ThreeVector(x_front_edge, y, z));
      updateMaxTheta(G4ThreeVector(x_back_edge, y, z));
    }

    // 2. Sample the vertices (center of caps)
    double x_front_vertex = lens.x - lens.tc / 2.0;
    double x_back_vertex  = lens.x + lens.tc / 2.0;
    updateMaxTheta(G4ThreeVector(x_front_vertex, 0, 0));
    updateMaxTheta(G4ThreeVector(x_back_vertex, 0, 0));

    // 3. Handle tangency to the spherical surfaces
    // xc1 is center of front surface sphere (facing -X)
    double xc1 = (lens.x - lens.tc / 2.0) + lens.R1;
    checkSphereTangency(sourcePos, axis, xc1, lens.R1, r, maxTheta);

    // xc2 is center of back surface sphere (facing +X)
    double xc2 = (lens.x + lens.tc / 2.0) - lens.R2;
    checkSphereTangency(sourcePos, axis, xc2, lens.R2, r, maxTheta);

    outMaxTheta = maxTheta;
  }

 private:
  static void checkSphereTangency(const G4ThreeVector& S, const G4ThreeVector& axis, double xc,
                                  double R, double r_lens, double& maxTheta) {
    if (R <= 0)
      return;

    // 2D reduction in the plane containing S and the X-axis
    double rs  = std::sqrt(S.y() * S.y() + S.z() * S.z());
    double dxs = xc - S.x();

    // Sphere center in 2D: (xc, 0)
    // Source in 2D: (S.x, rs)

    // Solve for tangent points (x, r) on the circle (x-xc)^2 + r^2 = R^2
    // dx^2 * (rs^2 + dxs^2) + dx * (2 * R^2 * dxs) + (R^4 - R^2 * rs^2) = 0

    double A     = rs * rs + dxs * dxs;
    double B     = 2.0 * R * R * dxs;
    double C_val = R * R * R * R - R * R * rs * rs;

    double delta = B * B - 4.0 * A * C_val;
    if (delta >= 0) {
      double sqrtDelta = std::sqrt(delta);
      double dx1       = (-B + sqrtDelta) / (2.0 * A);
      double dx2       = (-B - sqrtDelta) / (2.0 * A);

      for (double dx : {dx1, dx2}) {
        double x = xc + dx;
        // Check if this point is on the "visible" side of the sphere
        // and within the lens radius
        double r2 = R * R - dx * dx;
        if (r2 >= 0) {
          double r = std::sqrt(r2);
          // Check if this r is within the lens physical diameter
          if (r <= r_lens + 1e-6) {
            // In 3D, this tangent point corresponds to a circle on the sphere.
            // We need the point on this circle that maximizes the angle with 'axis'.
            // But wait, the 2D rs is already the max radial distance.
            // So we can just check the point (x, r) and (x, -r) in our 2D plane.

            // Convert back to 3D to use the provided 'axis'
            // The 2D plane is defined by the X-axis and the vector (0, Sy, Sz)
            G4ThreeVector rhoDir =
                (rs > 1e-9) ? G4ThreeVector(0, S.y(), S.z()).unit() : G4ThreeVector(0, 1, 0);

            auto check3D = [&](double px, double pr) {
              G4ThreeVector p   = G4ThreeVector(px, 0, 0) + pr * rhoDir;
              G4ThreeVector dir = (p - S).unit();
              double cosTheta   = dir.dot(axis);
              cosTheta          = std::max(-1.0, std::min(1.0, cosTheta));
              double theta      = std::acos(cosTheta);
              if (theta > maxTheta)
                maxTheta = theta;
            };

            check3D(x, r);
            check3D(x, -r);
          }
        }
      }
    }
  }
};

} // namespace riptide

#endif // RIPTIDE_IMPORTANCE_SAMPLING_HPP
