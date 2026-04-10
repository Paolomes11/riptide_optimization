
#include "importance_sampling.hpp"

#include <cmath>

namespace riptide {

void ImportanceSamplingHelper::CalculateGlobalCone(const std::vector<G4ThreeVector>& sourcePoints,
                                                   const LensParams& lens, G4ThreeVector& outAxis,
                                                   double& outMaxTheta) {
  if (sourcePoints.empty())
    return;

  G4ThreeVector avg(0, 0, 0);
  for (const auto& p : sourcePoints)
    avg += p;
  avg /= static_cast<double>(sourcePoints.size());

  G4ThreeVector lensCenter(lens.x, 0, 0);
  outAxis = (lensCenter - avg).unit();

  double maxTheta = 0.0;
  for (const auto& s : sourcePoints) {
    double t = 0.0;
    CalculateCone(s, lens, outAxis, t);
    if (t > maxTheta)
      maxTheta = t;
  }
  outMaxTheta = maxTheta;
}

void ImportanceSamplingHelper::CalculateCone(const G4ThreeVector& S, const LensParams& lens,
                                             const G4ThreeVector& axis, double& outMaxTheta) {
  double maxTheta = 0.0;
  auto update     = [&](const G4ThreeVector& p) {
    G4ThreeVector d = (p - S).unit();
    double c        = std::max(-1.0, std::min(1.0, d.dot(axis)));
    double theta    = std::acos(c);
    if (theta > maxTheta)
      maxTheta = theta;
  };

  const int nPhi = 48;
  double r       = lens.diameter * 0.5;

  double x_front_edge = lens.x - lens.te * 0.5;
  double x_back_edge  = lens.x + lens.te * 0.5;

  double x_front_vertex = lens.x - lens.tc * 0.5;
  double x_back_vertex  = lens.x + lens.tc * 0.5;

  for (int i = 0; i < nPhi; ++i) {
    double phi = (2.0 * M_PI * i) / nPhi;
    double y   = r * std::cos(phi);
    double z   = r * std::sin(phi);
    update(G4ThreeVector(x_front_edge, y, z));
    update(G4ThreeVector(x_back_edge, y, z));
    update(G4ThreeVector(x_front_vertex, y, z));
    update(G4ThreeVector(x_back_vertex, y, z));
  }

  update(G4ThreeVector(x_front_vertex, 0, 0));
  update(G4ThreeVector(x_back_vertex, 0, 0));

  if (lens.R1 > 0.0) {
    double xc1 = (lens.x - lens.tc * 0.5) + lens.R1;
    checkSphereTangency(S, axis, xc1, lens.R1, r, maxTheta);
  }
  if (lens.R2 > 0.0) {
    double xc2 = (lens.x + lens.tc * 0.5) - lens.R2;
    checkSphereTangency(S, axis, xc2, lens.R2, r, maxTheta);
  }

  double x_min = std::min(x_front_edge, x_front_vertex);
  double x_max = std::max(x_back_edge, x_back_vertex);
  checkCylinderTangency(S, axis, r, x_min, x_max, maxTheta);

  outMaxTheta = maxTheta;
}

void ImportanceSamplingHelper::checkSphereTangency(const G4ThreeVector& S,
                                                   const G4ThreeVector& axis, double xc, double R,
                                                   double r_lens, double& maxTheta) {
  if (R <= 0.0)
    return;

  double rs  = std::sqrt(S.y() * S.y() + S.z() * S.z());
  double dxs = xc - S.x();

  double A     = rs * rs + dxs * dxs;
  double B     = 2.0 * R * R * dxs;
  double C_val = R * R * R * R - R * R * rs * rs;

  double delta = B * B - 4.0 * A * C_val;
  if (delta < 0.0)
    return;

  double sqrtDelta = std::sqrt(delta);
  double dx1       = (-B + sqrtDelta) / (2.0 * A);
  double dx2       = (-B - sqrtDelta) / (2.0 * A);

  auto check_dx = [&](double dx) {
    double x  = xc + dx;
    double r2 = R * R - dx * dx;
    if (r2 < 0.0)
      return;
    double rl = std::sqrt(r2);
    if (rl > r_lens + 1e-6)
      return;

    G4ThreeVector rhoDir =
        (rs > 1e-9) ? G4ThreeVector(0, S.y(), S.z()).unit() : G4ThreeVector(0, 1, 0);

    auto upd = [&](double pr) {
      G4ThreeVector P = G4ThreeVector(x, 0, 0) + pr * rhoDir;
      G4ThreeVector d = (P - S).unit();
      double c        = std::max(-1.0, std::min(1.0, d.dot(axis)));
      double theta    = std::acos(c);
      if (theta > maxTheta)
        maxTheta = theta;
    };
    upd(rl);
    upd(-rl);
  };

  check_dx(dx1);
  check_dx(dx2);
}

void ImportanceSamplingHelper::checkCylinderTangency(const G4ThreeVector& S,
                                                     const G4ThreeVector& axis, double r,
                                                     double x_min, double x_max, double& maxTheta) {
  double rho = std::hypot(S.y(), S.z());
  if (rho <= r + 1e-9)
    return;

  double phi_c = std::atan2(S.z(), S.y());
  double delta = std::acos(r / rho);
  double phi1  = phi_c + delta;
  double phi2  = phi_c - delta;
  double y1    = r * std::cos(phi1);
  double z1    = r * std::sin(phi1);
  double y2    = r * std::cos(phi2);
  double z2    = r * std::sin(phi2);
  G4ThreeVector P1a(x_min, y1, z1), P1b(x_max, y1, z1);
  G4ThreeVector P2a(x_min, y2, z2), P2b(x_max, y2, z2);

  auto upd = [&](const G4ThreeVector& P) {
    G4ThreeVector d = (P - S).unit();
    double c        = std::max(-1.0, std::min(1.0, d.dot(axis)));
    double theta    = std::acos(c);
    if (theta > maxTheta)
      maxTheta = theta;
  };

  upd(P1a);
  upd(P1b);
  upd(P2a);
  upd(P2b);
}

} // namespace riptide
