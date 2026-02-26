#include "efficiency.hpp"
#include <fmt/core.h>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <TFile.h>
#include <TH2D.h>

constexpr double r1 = 38.6; // raggio lente 1 [mm]
constexpr double h1 = 12.5; // spessore lente 1 [mm]
constexpr double r2 = 30.9; // raggio lente 2 [mm]
constexpr double h2 = 16.3; // spessore lente 2 [mm]

constexpr double x_min = 1.0;
constexpr double x_max = 173.0;

constexpr double dx     = 3.0; // passo griglia [mm]
constexpr int n_photons = 10000;

const std::filesystem::path macro_path    = "run1.mac";
const std::filesystem::path output_file   = "output/output.root";
const std::filesystem::path root_map_file = "output/efficiency_map.root";

int main() {
  double best_efficiency = -1.0;
  double best_x1         = 0.0;
  double best_x2         = 0.0;

  double x1_min = x_min - r1 + h1;
  double x1_max = x_max - h2 - 1 - r1;

  // TH2D per mappa 2D: x1 su asse X, x2 su asse Y
  int nx_bins = static_cast<int>((x1_max - x1_min) / dx) + 1;
  int ny_bins = static_cast<int>((x_max - x_min) / dx) + 1; // max range approssimativo
  TH2D h_eff("h_eff", "Mappa efficienza; x1 [mm]; x2 [mm]; Efficienza", nx_bins, x1_min, x1_max,
             ny_bins, x_min, x_max);

  for (double x1 = x1_min; x1 <= x1_max; x1 += dx) {
    double x2_min = x1 + r1 + r2 + 1;
    double x2_max = x_max + r2 - h2;

    for (double x2 = x2_min; x2 <= x2_max; x2 += dx) {
      // Comando di simulazione
      std::string cmd =
          fmt::format("./build/Debug/simulate {:.3f} {:.3f} {}", x1, x2, macro_path.string());
      int ret = std::system(cmd.c_str());
      if (ret != 0) {
        std::cerr << fmt::format("Errore simulazione x1={:.1f}, x2={:.1f}\n", x1, x2);
        continue;
      }

      double efficiency = compute_geometric_efficiency(output_file.string(), n_photons);
      std::cout << fmt::format("x1={:.1f} mm, x2={:.1f} mm → efficienza={:.4f}\n", x1, x2,
                               efficiency);

      // Aggiorna best
      if (efficiency > best_efficiency) {
        best_efficiency = efficiency;
        best_x1         = x1;
        best_x2         = x2;
      }

      // Salva nella mappa
      h_eff.Fill(x1, x2, efficiency);
    }
  }

  std::cout << fmt::format(
      "\nConfigurazione ottimale:\n  x1={:.2f} mm\n  x2={:.2f} mm\n  efficienza={:.4f}\n", best_x1,
      best_x2, best_efficiency);

  // Salva la mappa in un file ROOT
  TFile f(root_map_file.string().c_str(), "RECREATE");
  h_eff.Write();
  f.Close();

  std::cout << "Mappa di efficienza salvata in: " << root_map_file << "\n";

  return 0;
}