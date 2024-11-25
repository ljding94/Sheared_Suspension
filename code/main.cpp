// Copyright[2024] [Lijie Ding]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>
#include <filesystem>
#include "suspension.h"

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
    double Rmu = 1*R0;
    double sqrtD0 = 1*R0; // for simplicity: dx, dy = sqrt(2D0)*rand_norm(gen)
    // size dependence, Einstein: D ~ 1/a, a is particle diameter

    // precision run with specified parameters
    if (argc == 7)
    {
        int n = std::atoi(argv[1]); // number density
        double Rmu = std::atof(argv[2])*R0;
        double sigma = std::atof(argv[3]);
        double sqrtD = std::atof(argv[4])*R0;
        // double theta = std::atof(argv[3]);
        // double Sx = std::atof(argv[4]);
        // double phi = std::atof(argv[5]);
        double gxy = std::atof(argv[5])*R0;
        std::string folder = std::string(argv[6]);

        suspension gas_2d(R0, Rmu, n, sigma, sqrtD, gxy, false);

        std::string finfo = "n" + std::string(argv[1]) + "_Rmu" + std::string(argv[2]) + "_sigma" + std::string(argv[3]) + "_sqrtD" + std::string(argv[4]) + "_gxy" + std::string(argv[5]);

        int number_of_config = 2000;
        int bnum_r = 100;
        int bnum_phi = 101;

        gas_2d.run_simulation(number_of_config, bnum_r, bnum_phi, folder, finfo);
    }

    // random run
    if (argc == 3)
    {
        int run_num = std::atoi(argv[1]);
        std::string folder = std::string(argv[2]);

        // number_of_config 2000, bin_num 101, n 200 takes 468s to run

        suspension gas_2d(R0, 1.0*R0, 100, 0.0, 1.0, 0.0, true);
        std::string finfo = "random_run" + std::to_string(run_num);

        int number_of_config = 4000;
        int bnum_r = 100;
        int bnum_phi = 101;
        gas_2d.run_simulation(number_of_config, bnum_r, bnum_phi, folder, finfo);
    }

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << time_elapsed << " seconds" << std::endl;

    return 0;
}