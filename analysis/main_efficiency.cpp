#include "efficiency.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Uso: " << argv[0] << " <ROOT file> <numero fotoni sparati>\n";
        return 1;
    }

    std::string root_file = argv[1];
    int n_photons_shot = std::atoi(argv[2]);
    if (n_photons_shot <= 0)
    {
        std::cerr << "Errore: numero fotoni sparati deve essere > 0\n";
        return 1;
    }

    double efficiency = compute_geometric_efficiency(root_file, n_photons_shot);
    std::cout << "Efficienza geometrica calcolata: " << efficiency << std::endl;

    return 0;
}