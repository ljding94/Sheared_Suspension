#ifndef _IDEAL_GAS_H
#define _IDEAL_GAS_H
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

// observable
struct observable
{
    std::vector<double> qD{};                   // qB for Sq, D is the bead diameter, which is 1 in our unit
    std::vector<std::vector<double>> Sq2D{};    // structure factor
    std::vector<std::vector<double>> Sq2D_af{}; // structure factor after affinr transformation
    std::vector<std::vector<double>> SqSq2D{};  // Sq-Sq correlation
};
struct bead
{
    std::vector<double> r; // position of the bead
    double f;              // scattering weight, for log-normal ln(f) ~ N(0,sigma), PDI = e^{sigma^2}, sigma = sqrt(ln(PDI))
};

class ideal_gas
{
public:
    int n; // number of beads or number density, since use box size L=1

    double sigma;  // std of log-normal distribution of scattering weight
    double theta;  // first rotation angle
    double Sx, Sy; // streching factor
    double phi;    // second rotation angle

    double gamma_xx, gamma_xy, gamma_yx, gamma_yy; // affine matrix

    std::vector<bead> all_beads; // the gas sysem

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_real_distribution<> rand_uni; // uniform distribution
    std::normal_distribution<> rand_norm;     // normal distribution

    // initialization
    ideal_gas(int n_, double sigma_, double theta_, double Sx_, double phi_, bool random_param = false);

    int generate_gas();

    // some observable measurement
    observable measure_observable(int bin_num);

    void affine_transform();
    std::vector<std::vector<double>> calc_structure_factor_2d(std::vector<double> qD);
    std::vector<std::vector<double>> calc_structure_factor_2d_af(std::vector<double> qD);
    std::vector<std::vector<double>> calc_structure_factor_correlation(std::vector<std::vector<double>> Sq1, std::vector<std::vector<double>> Sq2, std::vector<double> qD);

    // experiment
    void save_gas_config_to_file(std::string filename);
    void save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble);
    void run_simulation(int N_config, int bin_num, std::string folder, std::string finfo);
};

#endif