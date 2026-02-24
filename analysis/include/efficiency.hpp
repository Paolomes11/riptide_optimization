#ifndef EFFICIENCY_HPP
#define EFFICIENCY_HPP

#include <string>
#include <vector>

// Calcola l'efficienza geometrica
// input: nome file ROOT, numero di fotoni sparati
// output: efficienza [0,1]
double compute_geometric_efficiency(const std::string&, int);

#endif // EFFICIENCY_HPP