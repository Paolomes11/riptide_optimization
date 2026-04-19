#include <nlohmann/json.hpp>

#include <TCanvas.h>
#include <TColor.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct CliConfig {
  std::string tsv_file    = "output/dof_map.tsv";
  std::string config_file = "config/config.json";
  std::string output_dir  = "output/dof_analysis";
  std::optional<double> m_target;
};

static CliConfig parse_args(int argc, char** argv) {
  CliConfig cfg;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    auto next       = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "Argomento mancante dopo " << arg << "\n";
        std::exit(1);
      }
      return argv[++i];
    };

    if (arg == "--tsv")
      cfg.tsv_file = next();
    else if (arg == "--config" || arg == "-c")
      cfg.config_file = next();
    else if (arg == "--output" || arg == "-o")
      cfg.output_dir = next();
    else if (arg == "--m-target")
      cfg.m_target = std::stod(next());
    else if (arg == "--help" || arg == "-h") {
      std::cout << "Uso:\n"
                << "  magnification_map --tsv output/dof_map.tsv [--m-target 0.1333]\n"
                << "                  [--config config/config.json] [--output dir/]\n";
      std::exit(0);
    } else {
      std::cerr << "Opzione sconosciuta: " << arg << "\n";
      std::exit(1);
    }
  }
  return cfg;
}

static void set_root_style() {
  gStyle->Reset();
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleFont(42, "");
  gStyle->SetStatFont(42);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNumberContours(255);
}

static std::vector<std::string> split_tab(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  cur.reserve(s.size());
  for (char c : s) {
    if (c == '\t') {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

static std::string trim_cr(const std::string& s) {
  if (!s.empty() && s.back() == '\r') {
    return s.substr(0, s.size() - 1);
  }
  return s;
}

struct Row {
  double x1 = 0.0;
  double x2 = 0.0;
  double M  = 0.0;
};

static double min_positive_step(const std::vector<double>& v_sorted) {
  double dmin = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < v_sorted.size(); ++i) {
    double d = v_sorted[i] - v_sorted[i - 1];
    if (d > 0.0 && d < dmin) {
      dmin = d;
    }
  }
  if (!std::isfinite(dmin)) {
    return 1.0;
  }
  return dmin;
}

static void set_diverging_palette_red_white_green() {
  const int nRGBs     = 3;
  double stops[nRGBs] = {0.0, 0.5, 1.0};
  double red[nRGBs]   = {1.0, 1.0, 0.0};
  double green[nRGBs] = {0.0, 1.0, 1.0};
  double blue[nRGBs]  = {0.0, 1.0, 0.0};
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
}

int main(int argc, char** argv) {
  using json = nlohmann::json;

  CliConfig cli = parse_args(argc, argv);
  set_root_style();
  std::filesystem::create_directories(cli.output_dir);

  std::ifstream jf(cli.config_file);
  if (!jf.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.config_file << "\n";
    return 1;
  }
  json config;
  jf >> config;

  double m_target = cli.m_target.value_or(config.value("m_target", 1.0 / 7.5));

  std::ifstream f(cli.tsv_file);
  if (!f.is_open()) {
    std::cerr << "Errore: impossibile aprire " << cli.tsv_file << "\n";
    return 1;
  }

  std::string header_line;
  if (!std::getline(f, header_line)) {
    std::cerr << "Errore: TSV vuoto\n";
    return 1;
  }
  auto header = split_tab(trim_cr(header_line));
  std::unordered_map<std::string, size_t> col;
  col.reserve(header.size());
  for (size_t i = 0; i < header.size(); ++i) {
    col[header[i]] = i;
  }

  auto get_idx = [&](const std::string& name) -> size_t {
    auto it = col.find(name);
    if (it == col.end()) {
      std::cerr << "Errore: colonna TSV mancante: " << name << "\n";
      std::exit(1);
    }
    return it->second;
  };

  size_t ix1 = get_idx("x1");
  size_t ix2 = get_idx("x2");
  size_t iM  = get_idx("M");

  std::vector<Row> rows;
  rows.reserve(4096);

  std::vector<double> x1_vals;
  std::vector<double> x2_vals;
  x1_vals.reserve(4096);
  x2_vals.reserve(4096);

  std::string line;
  while (std::getline(f, line)) {
    line = trim_cr(line);
    if (line.empty()) {
      continue;
    }
    auto fields = split_tab(line);
    size_t need = std::max({ix1, ix2, iM}) + 1;
    if (fields.size() < need) {
      continue;
    }
    Row r;
    r.x1 = std::stod(fields[ix1]);
    r.x2 = std::stod(fields[ix2]);
    r.M  = std::stod(fields[iM]);
    rows.push_back(r);
    x1_vals.push_back(r.x1);
    x2_vals.push_back(r.x2);
  }

  if (rows.empty()) {
    std::cerr << "Errore: nessuna riga dati nel TSV\n";
    return 1;
  }

  std::sort(x1_vals.begin(), x1_vals.end());
  x1_vals.erase(std::unique(x1_vals.begin(), x1_vals.end()), x1_vals.end());
  std::sort(x2_vals.begin(), x2_vals.end());
  x2_vals.erase(std::unique(x2_vals.begin(), x2_vals.end()), x2_vals.end());

  if (x1_vals.size() < 1 || x2_vals.size() < 1) {
    std::cerr << "Errore: griglia (x1,x2) non valida\n";
    return 1;
  }

  double dx1 = (x1_vals.size() >= 2) ? min_positive_step(x1_vals) : 1.0;
  double dx2 = (x2_vals.size() >= 2) ? min_positive_step(x2_vals) : 1.0;

  int n_bins_x1 = static_cast<int>(x1_vals.size());
  int n_bins_x2 = static_cast<int>(x2_vals.size());
  double x1_min = x1_vals.front();
  double x1_max = x1_vals.back();
  double x2_min = x2_vals.front();
  double x2_max = x2_vals.back();

  TH2D h_M_error("h_M_error", ";x1 [mm];x2 [mm];M - M_{target}", n_bins_x1, x1_min - dx1 / 2.0,
                 x1_max + dx1 / 2.0, n_bins_x2, x2_min - dx2 / 2.0, x2_max + dx2 / 2.0);
  TH2D h_M_abs_err("h_M_abs_err", ";x1 [mm];x2 [mm];|M - M_{target}|", n_bins_x1,
                   x1_min - dx1 / 2.0, x1_max + dx1 / 2.0, n_bins_x2, x2_min - dx2 / 2.0,
                   x2_max + dx2 / 2.0);

  double best_abs = std::numeric_limits<double>::infinity();
  std::pair<double, double> best_xy{0.0, 0.0};

  double max_abs_err = 0.0;
  for (const auto& r : rows) {
    double err     = r.M - m_target;
    double abs_err = std::abs(err);
    int bx         = h_M_error.GetXaxis()->FindBin(r.x1);
    int by         = h_M_error.GetYaxis()->FindBin(r.x2);
    h_M_error.SetBinContent(bx, by, err);
    h_M_abs_err.SetBinContent(bx, by, abs_err);
    max_abs_err = std::max(max_abs_err, abs_err);
    if (abs_err < best_abs) {
      best_abs = abs_err;
      best_xy  = {r.x1, r.x2};
    }
  }

  if (!(max_abs_err > 0.0)) {
    max_abs_err = 1.0;
  }
  h_M_error.SetMinimum(-max_abs_err);
  h_M_error.SetMaximum(+max_abs_err);

  auto save_diverging = [&](TH2D& h, const std::string& name) {
    set_diverging_palette_red_white_green();
    TCanvas c(("c_" + name).c_str(), name.c_str(), 1100, 900);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.08);
    h.Draw("COLZ");
    TMarker m(best_xy.first, best_xy.second, 29);
    m.SetMarkerSize(2.0);
    m.SetMarkerColor(kBlack);
    m.Draw("same");
    std::string out = (std::filesystem::path(cli.output_dir) / (name + ".png")).string();
    c.SaveAs(out.c_str());
  };

  auto save_viridis = [&](TH2D& h, const std::string& name) {
    gStyle->SetPalette(kViridis);
    TCanvas c(("c_" + name).c_str(), name.c_str(), 1100, 900);
    c.SetLeftMargin(0.16);
    c.SetBottomMargin(0.14);
    c.SetRightMargin(0.16);
    c.SetTopMargin(0.08);
    h.Draw("COLZ");
    TMarker m(best_xy.first, best_xy.second, 29);
    m.SetMarkerSize(2.0);
    m.SetMarkerColor(kBlack);
    m.Draw("same");
    std::string out = (std::filesystem::path(cli.output_dir) / (name + ".png")).string();
    c.SaveAs(out.c_str());
  };

  save_diverging(h_M_error, "magnification_M_error_map");
  save_viridis(h_M_abs_err, "magnification_M_abs_error_map");

  return 0;
}
