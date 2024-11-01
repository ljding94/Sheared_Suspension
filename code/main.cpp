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

    // precision run with specified parameters
    if (argc == 7)
    {
        int n = std::atoi(argv[1]); // number density
        double sigma = std::atof(argv[2]);
        double theta = std::atof(argv[3]);
        double Sx = std::atof(argv[4]);
        double phi = std::atof(argv[5]);
        std::string folder = std::string(argv[6]);

        ideal_gas gas_2d(n, sigma, theta, Sx, phi);
        // std::string finfo = "L" + std::string(argv[1]) + "_N" + std::string(argv[2]) + "_theta" + std::string(argv[3]) + "_Sx" + std::string(argv[4]) + "_Sy" + std::string(argv[5]) + "_phi" + std::string(argv[6]);
        std::string finfo = "n" + std::string(argv[1]) + "_sigma" + std::string(argv[2]) + "_theta" + std::string(argv[3]) + "_Sx" + std::string(argv[4]) + "_phi" + std::string(argv[5]);

        int number_of_config = 2000;
        int bin_num = 51;

        gas_2d.run_simulation(number_of_config, bin_num, folder, finfo);
    }

    // random run
    if (argc == 4)
    {
        int n = std::atoi(argv[1]);
        int run_num = std::atoi(argv[2]);
        std::string folder = std::string(argv[3]);

        // number_of_config 2000, bin_num 101, n 200 takes 468s to run

        int number_of_config = 2000;
        int bin_num = 101;

        ideal_gas gas_2d(n, 0, 0, 1, 0, true); // random theta,Sx,Sy,phi
        std::string finfo = std::string(argv[1]) + "_random_run" + std::to_string(run_num);

        gas_2d.run_simulation(number_of_config, bin_num, folder, finfo);
    }

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << time_elapsed << " seconds" << std::endl;

    return 0;
}