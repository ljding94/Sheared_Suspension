// Copyright[2024] [Lijie Ding]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>
#include <filesystem>
#include "ideal_gas.h"

int main(int argc, char const *argv[])
{
    std::clock_t c_start = std::clock();
    std::string folder;
    double beta = 1;

    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm *timeinfo = std::localtime(&now);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);
    std::string today(buffer);
    std::cout << today << std::endl;

    std::cout << "running with input argc:" << argc << "\n";
    for (int i = 0; i < argc; i++)
    {
        std::cout << argv[i] << " ";
    }

    double R0 = 0.00125;
    double sqrtD0 = 1*R0; // for simplicity: dx, dy = sqrt(2D0)*rand_norm(gen)
    // size dependence, Einstein: D ~ 1/a, a is particle diameter

    // precision run with specified parameters
    if (argc == 8)
    {
        int N = std::atoi(argv[1]); // number density
        double sigma = std::atof(argv[2]);
        double sqrtD0 = 1*R0;
        // double theta = std::atof(argv[3]);
        // double Sx = std::atof(argv[4]);
        // double phi = std::atof(argv[5]);

        double gxx = std::atof(argv[3]);
        double gxy = std::atof(argv[4]);
        double gyx = std::atof(argv[5]);
        double gyy = std::atof(argv[6]);
        std::string folder = std::string(argv[7]);

        ideal_gas gas_2d(N, R0, sqrtD0, sigma, gxx, gxy, gyx, gyy);
        // ideal_gas gas_2d(N, sigma, theta, Sx, phi);
        //  std::string finfo = "L" + std::string(argv[1]) + "_N" + std::string(argv[2]) + "_theta" + std::string(argv[3]) + "_Sx" + std::string(argv[4]) + "_Sy" + std::string(argv[5]) + "_phi" + std::string(argv[6]);
        // std::string finfo = "N" + std::string(argv[1]) + "_sigma" + std::string(argv[2]) + "_theta" + std::string(argv[3]) + "_Sx" + std::string(argv[4]) + "_phi" + std::string(argv[5]);
        std::string finfo = "N" + std::string(argv[1]) + "_sigma" + std::string(argv[2]) + "_gxx" + std::string(argv[3]) + "_gxy" + std::string(argv[4]) + "_gyx" + std::string(argv[5]) + "_gyy" + std::string(argv[6]);

        int number_of_config = 2000;
        int bnum_r = 100;
        int bnum_phi = 101;

        gas_2d.run_simulation(number_of_config, bnum_r, bnum_phi, folder, finfo);
    }

    // random run
    if (argc == 4)
    {
        int N = std::atoi(argv[1]);
        int run_num = std::atoi(argv[2]);
        std::string folder = std::string(argv[3]);

        // number_of_config 2000, bin_num 101, n 200 takes 468s to run

        ideal_gas gas_2d(N, R0, sqrtD0, 0, 1, 0, 0, 1, true); // random theta,Sx,Sy,phi
        std::string finfo = "N" + std::string(argv[1]) + "_random_run" + std::to_string(run_num);

        int number_of_config = 20000;
        int bnum_r = 50;
        int bnum_phi = 51;
        gas_2d.run_simulation(number_of_config, bnum_r, bnum_phi, folder, finfo);
    }

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << time_elapsed << " seconds" << std::endl;

    return 0;
}