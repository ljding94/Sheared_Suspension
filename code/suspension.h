#ifndef _SUSPENSION_H
#define _SUSPENSION_H
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

// observable
struct observable
{
    std::vector<double> qr{};
    std::vector<double> qphi{};

    std::vector<std::vector<double>> Iq2D{};    // scattering function
    std::vector<std::vector<double>> Iq2D_af{}; // scartering function after affinr transformation
    std::vector<std::vector<double>> IqIq_af{};  // scattering function correlation
    //TODO: get back to this and continue
};
struct bead
{
    std::vector<double> r; // position of the bead
    double R;              // 0.025, on average radius of the bead
    double V;              // volume of the bead
    std::vector<double> fq; // form factor of the bead
    // for sphereical particle:
    // fq = 3 * (sin(qR) - qR cos(qR)) / (qR)^3 when q>0,
    // fq = 1 when q=0
    // form factor: http://gisaxs.com/index.php/Form_Factor:Sphere
    // note: Fq = |fq|^2
};

class suspension
{
public:
    double R0;    // useful length scale
    double Rmu;    // average radius of the bead
    int n; // number of beads or number density, since use box size L=1
    double sigma;  // std of log-normal distribution of scattering weight
    double sqrtD;    // average diffusion coefficient
    double gxy; // shear

    std::vector<bead> all_beads; // the gas sysem

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_real_distribution<> rand_uni; // uniform distribution
    std::normal_distribution<> rand_norm;      // normal distribution

    // initialization
    //suspension(int n_, double sigma_, double theta_, double Sx_, double phi_, bool random_param = false);
    suspension(double R0_, double Rmu_, int n_, double sigma_, double sqrtD_, double gxy_, bool random_param = false);

    int generate_gas();

    // some observable measurement
    observable measure_observable(int bnum_r, int bnum_phi);

    void affine_transform();
    void Brownian_transform();
    std::vector<double> calc_structure_sphere_form_factor(std::vector<double> qr, double R);
    std::vector<std::vector<double>> calc_structure_factor_2d(std::vector<double> qr,std::vector<double> qphi);
    std::vector<std::vector<double>> calc_structure_factor_2d_af(std::vector<double> qr,std::vector<double> qphi);

    std::vector<std::vector<double>> calc_structure_factor_correlation(std::vector<std::vector<double>> Sq1, std::vector<std::vector<double>> Sq2, std::vector<double> qD);


    std::vector<std::vector<double>> element_wise_dot_2d(std::vector<std::vector<double>> M1, std::vector<std::vector<double>> M2);

    // experiment
    void save_gas_config_to_file(std::string filename);
    void save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble);
    void run_simulation(int N_config, int bnum_r, int bnum_phi, std::string folder, std::string finfo);
};

#endif