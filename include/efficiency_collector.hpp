#ifndef RIPTIDE_EFFICIENCY_COLLECTOR_HPP
#define RIPTIDE_EFFICIENCY_COLLECTOR_HPP

#include <set>

namespace riptide {

class EfficiencyCollector {
  std::set<int> m_events_with_hits;

 public:
  // registra un evento che ha prodotto almeno un hit
  void recordEvent(int event_id);

  // registra un evento che ha prodotto almeno un hit (metodo corretto)
  void recordEventWithHit(int event_id);

  // calcola l'efficienza geometrica
  double computeEfficiency(int n_photons_shot) const;

  // resetta i dati tra un run e l'altro
  void reset();

  // Accesso globale al collector (per il sensitive detector)
  static EfficiencyCollector* GetInstance();
  static void SetInstance(EfficiencyCollector* collector);

 private:
  static EfficiencyCollector* s_instance;
};

} // namespace riptide

#endif // RIPTIDE_EFFICIENCY_COLLECTOR_HPP