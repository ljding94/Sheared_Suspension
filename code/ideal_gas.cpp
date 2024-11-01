#include "ideal_gas.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>

// initialization
ideal_gas::ideal_gas(int n_, double sigma_, double theta_, double Sx_, double phi_, bool random_param)
{

    // set random number generators
    std::random_device rd;
    std::mt19937 gen_set(rd());
    std::uniform_real_distribution<> rand_uni_set(0, 1);
    std::normal_distribution<> rand_norm_set(0, 1); //

    gen = gen_set;
    rand_uni = rand_uni_set;
    rand_norm = rand_norm_set;

    n = n_;

    // generate the gas system
    if (random_param)
    {
        sigma = 1.0 * rand_uni(gen);
        theta = 0.1*M_PI * rand_uni(gen);
        Sx = 0.9 + 0.2 * rand_uni(gen);
        phi = 0.1*M_PI * rand_uni(gen);
    }
    else
    {
        sigma = sigma_;
        theta = theta_;
        Sx = Sx_;
        phi = phi_;
    }
    Sy = 1.0 / Sx;
    std::cout << "theta: " << theta << ", Sx: " << Sx << ", Sy: " << Sy << ", phi: " << phi << std::endl;

    // calc the affine matrix
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    gamma_xx = Sx * cos_theta * cos_phi + Sy * sin_theta * sin_phi;
    gamma_xy = Sx * sin_theta * cos_phi - Sy * cos_theta * sin_phi;
    gamma_yx = Sx * cos_theta * sin_phi - Sy * sin_theta * cos_phi;
    gamma_yy = Sx * sin_theta * sin_phi + Sy * cos_theta * cos_phi;
}

int ideal_gas::generate_gas()
{
    // for simplicity, let's ignore hard disk interaction and just direct sample random positions
    // generate the gas system
    // TODO: generate gas in a larger box, and always calculate Sq for particles in the smaller sqaure box, for consistency with experiments
    all_beads.resize(4 * n); // generate beads in 2x2 squre, measure Sq for 1x1 square
    double r_theta;
    double r_r;
    double L = 2.0; // box size
    for (int i = 0; i < all_beads.size(); i++)
    {
        all_beads[i].f = std::exp(sigma * rand_norm(gen)); // log-normal distribution, PDI=e^{sigma^2}
        all_beads[i].r = {(rand_uni(gen) - 0.5) * L, (rand_uni(gen) - 0.5) * L};
        //  let's try circular shape
        // r_r = L * std::sqrt(rand_uni(gen));
        // r_theta = 2 * M_PI * rand_uni(gen);
        // all_beads[i].r = {r_r * std::cos(r_theta), r_r * std::sin(r_theta)};
    }
    return 1;
}

observable ideal_gas::measure_observable(int bin_num)
{
    // measure the observable
    observable obs;

    // measure the structure factor
    double qD_i = -25.0 * M_PI;
    double dqD = 50.0 * M_PI / (bin_num - 1);
    obs.qD.resize(bin_num);
    for (int k = 0; k < bin_num; k++)
    {
        obs.qD[k] = qD_i + dqD * k;
    }
    obs.Sq2D = calc_structure_factor_2d(obs.qD);

    obs.Sq2D_af = calc_structure_factor_2d_af(obs.qD);

    // measure the Sq-Sq correlation
    std::vector<std::vector<double>> SqSq(bin_num, std::vector<double>(bin_num, 0.0));
    std::vector<std::vector<double>> SqafSqaf(bin_num, std::vector<double>(bin_num, 0.0));
    std::vector<std::vector<double>> SqSqaf(bin_num, std::vector<double>(bin_num, 0.0));
    std::vector<std::vector<double>> SqSqaf_norm(bin_num, std::vector<double>(bin_num, 0.0));

    SqSq = calc_structure_factor_correlation(obs.Sq2D, obs.Sq2D, obs.qD);
    SqafSqaf = calc_structure_factor_correlation(obs.Sq2D_af, obs.Sq2D_af, obs.qD);
    SqSqaf = calc_structure_factor_correlation(obs.Sq2D, obs.Sq2D_af, obs.qD);

    for (int kx = 0; kx < bin_num; kx++)
    {
        for (int ky = 0; ky < bin_num; ky++)
        {
            SqSqaf_norm[kx][ky] = SqSqaf[kx][ky] / std::sqrt(SqSq[kx][ky] * SqafSqaf[kx][ky]);
        }
    }
    obs.SqSq2D = SqSqaf_norm;

    return obs;
}

std::vector<std::vector<double>> ideal_gas::calc_structure_factor_2d(std::vector<double> qD)
{
    // measrue 2d structure factor
    int bin_num = qD.size();
    int bin0 = (bin_num - 1) / 2;
    if (bin_num % 2 == 0)
    {
        std::cout << "Error: bin_num should be odd" << std::endl;
    }
    if (std::abs(qD[bin0]) > 1e-6)
    {
        std::cout << "Error: qD[bin0] should be 0" << std::endl;
        std::cout << "qD[bin0]" << qD[bin0] << std::endl;
    }

    std::vector<double> r{0, 0};
    double qx, qy, Sq_Re_buff, Sq_Im_buff;
    std::vector<std::vector<double>> Sq(bin_num, std::vector<double>(bin_num, 0.0));
    std::vector<std::vector<double>> Sq_Re(bin_num, std::vector<double>(bin_num, 0.0));
    std::vector<std::vector<double>> Sq_Im(bin_num, std::vector<double>(bin_num, 0.0));

    std::vector<bead> scatter_beads;
    scatter_beads.resize(0);

    // add all in beam square range beads
    for (int i = 0; i < all_beads.size(); i++)
    {
        if (all_beads[i].r[0] < 0.5 && all_beads[i].r[0] > -0.5 && all_beads[i].r[1] < 0.5 && all_beads[i].r[1] > -0.5)
        {
            scatter_beads.push_back(all_beads[i]);
        }
    }
    int N_scatter = scatter_beads.size();

    // calculate S_q
    for (int kx = 0; kx < bin_num; kx++)
    {
        qx = qD[kx];
        for (int ky = bin0; ky < bin_num; ky++)
        {
            for (int i = 0; i < scatter_beads.size() - 1; i++)
            {
                r[0] = scatter_beads[i].r[0];
                r[1] = scatter_beads[i].r[1];
                qy = qD[ky];
                Sq_Re_buff = 1.0 * scatter_beads[i].f / N_scatter * std::cos(qx * r[0] + qy * r[1]);
                Sq_Re[kx][ky] += Sq_Re_buff;
                Sq_Im_buff = 1.0 * scatter_beads[i].f / N_scatter * std::sin(qx * r[0] + qy * r[1]);
                Sq_Im[kx][ky] += Sq_Im_buff;
                if (ky != bin0)
                {
                    Sq_Re[2 * bin0 - kx][2 * bin0 - ky] += Sq_Re_buff;
                    Sq_Im[2 * bin0 - kx][2 * bin0 - ky] += Sq_Im_buff;
                }
            }
            Sq[kx][ky] = Sq_Re[kx][ky] * Sq_Re[kx][ky] + Sq_Im[kx][ky] * Sq_Im[kx][ky];
            if (ky != bin0)
            {
                Sq[2 * bin0 - kx][2 * bin0 - ky] = Sq[kx][ky];
            }
        }
    }
    return Sq;
}

void ideal_gas::affine_transform()
{
    // affine transform the gas system
    double x, y;
    for (int i = 0; i < all_beads.size(); i++)
    {
        x = all_beads[i].r[0];
        y = all_beads[i].r[1];
        all_beads[i].r[0] = gamma_xx * x + gamma_xy * y;
        all_beads[i].r[1] = gamma_yx * x + gamma_yy * y;
    }
}

std::vector<std::vector<double>> ideal_gas::calc_structure_factor_2d_af(std::vector<double> qD)
{
    // measrue 2d structure factor
    int bin_num = qD.size();
    std::vector<std::vector<double>> Sq; // initialize with 1 due to self overlaping term (see S.H. Chen 1986 eq 18)
    std::vector<bead> all_beads_pre_af = all_beads;
    affine_transform();
    Sq = calc_structure_factor_2d(qD);
    all_beads = all_beads_pre_af;
    return Sq;
}

std::vector<std::vector<double>> ideal_gas::calc_structure_factor_correlation(std::vector<std::vector<double>> Sq1, std::vector<std::vector<double>> Sq2, std::vector<double> qD)
{
    // TODO: normalize by autocorrelations

    // for simplicity let qD = deldD
    int bin_num = qD.size();
    int bin0 = (bin_num - 1) / 2;

    std::vector<double> r{0, 0};
    double delqx, delqy, SqSq_buff;
    std::vector<std::vector<double>> SqSq(bin_num, std::vector<double>(bin_num, 0.0));
    double SqSq0; // del_qx and del_qy = 0, for normalization
    double qD_max = qD[qD.size() - 1];
    int dqx_count = 0;
    int dqy_count = 0;
    int xi, xf, yi, yf;

    for (int dkx = 0; dkx < bin_num; dkx++)
    {
        for (int dky = 0; dky < bin_num; dky++)
        {
            // for each del_q, calculate the SqSq
            SqSq_buff = 0;
            // loop over all qx, qy;
            xi = std::max(0, dkx - bin0);
            xf = std::min(bin0 + dkx + 1, bin_num);
            yi = std::max(0, dky - bin0);
            yf = std::min(bin0 + dky + 1, bin_num);
            for (int i = xi; i < xf; i++)
            {
                for (int j = yi; j < yf; j++)
                {
                    SqSq[dkx][dky] += Sq1[i][j] * Sq2[i - dkx + bin0][j - dky + bin0];
                }
            }
            dqx_count = std::min(bin0 + dkx + 1, bin_num + bin0 - dkx);
            dqy_count = std::min(bin0 + dky + 1, bin_num + bin0 - dky);
            SqSq[dkx][dky] /= dqx_count * dqy_count;
        }
    }
    SqSq0 = SqSq[bin0][bin0];
    for (int dkx = 0; dkx < bin_num; dkx++)
    {
        for (int dky = 0; dky < bin_num; dky++)
        {
            SqSq[dkx][dky] /= SqSq0;
        }
    }
    return SqSq;
}

void ideal_gas::save_gas_config_to_file(std::string filename)
{
    // save the gas system to file
    std::ofstream f(filename);
    if (f.is_open())
    {
        std::vector<bead> all_beads_pre_af = all_beads;
        affine_transform();
        f << "f,x,y,x_af,y_af\n";
        for (int i = 0; i < all_beads.size(); i++)
        {
            f << all_beads_pre_af[i].f << "," << all_beads_pre_af[i].r[0] << "," << all_beads_pre_af[i].r[1] << "," << all_beads[i].r[0] << "," << all_beads[i].r[1] << "\n";
        }
        f.close();
    }
    else
    {
        std::cout << "Error: cannot open file" << filename << std::endl;
    }
}

void ideal_gas::save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble)
{
    // save the observable to file
    std::ofstream f(filename);
    if (f.is_open())
    {
        int number_of_config = obs_ensemble.size();
        // find stats
        // calculate average and standard deviation of E
        std::vector<std::vector<double>> avg_Sq2D(obs_ensemble[0].Sq2D.size(), std::vector<double>(obs_ensemble[0].Sq2D[0].size(), 0.0));
        std::vector<std::vector<double>> avg_Sq2D_af(obs_ensemble[0].Sq2D.size(), std::vector<double>(obs_ensemble[0].Sq2D[0].size(), 0.0));
        std::vector<std::vector<double>> avg_SqSq2D(obs_ensemble[0].SqSq2D.size(), std::vector<double>(obs_ensemble[0].SqSq2D[0].size(), 0.0));

        double M = obs_ensemble.size();
        double sqrt_M = std::sqrt(M);

        // get the average
        for (int i = 0; i < obs_ensemble.size(); i++)
        {
            for (int kx = 0; kx < obs_ensemble[i].Sq2D.size(); kx++)
            {
                for (int ky = 0; ky < obs_ensemble[i].Sq2D[kx].size(); ky++)
                {
                    avg_Sq2D[kx][ky] += obs_ensemble[i].Sq2D[kx][ky];
                    avg_Sq2D_af[kx][ky] += obs_ensemble[i].Sq2D_af[kx][ky];
                    avg_SqSq2D[kx][ky] += obs_ensemble[i].SqSq2D[kx][ky];
                }
            }
        }

        for (int kx = 0; kx < obs_ensemble[0].Sq2D.size(); kx++)
        {
            for (int ky = 0; ky < obs_ensemble[0].Sq2D[kx].size(); ky++)
            {
                avg_Sq2D[kx][ky] /= M;
                avg_Sq2D_af[kx][ky] /= M;
                avg_SqSq2D[kx][ky] /= M;
            }
        }

        f << "label,n,theta,Sx,Sy,phi,gamma_xx,gamma_xy,gamma_yz,gamma_yy,Sq2D/SqSq2D\n";
        // write parameters and stats to the file
        f << "mean," << n << "," << theta << "," << Sx << "," << Sy << "," << phi
         << "," << gamma_xx << "," << gamma_xy << "," << gamma_yx << "," << gamma_yy;
        for (int kx = 0; kx < obs_ensemble[0].Sq2D.size(); kx++)
        {
            f << ",NA";
        }
        f << "\n qD,NA,NA,NA,NA,NA,NA,NA,NA,NA";
        for (int j = 0; j < obs_ensemble[0].qD.size(); j++)
        {
            f << "," << obs_ensemble[0].qD[j];
        }

        for (int kx = 0; kx < obs_ensemble[0].Sq2D.size(); kx++)
        {
            f << "\nSq2D,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int ky = 0; ky < obs_ensemble[0].Sq2D[kx].size(); ky++)
            {
                f << "," << avg_Sq2D[kx][ky];
            }
        }

        for (int kx = 0; kx < obs_ensemble[0].Sq2D.size(); kx++)
        {
            f << "\nSq2D_af,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int ky = 0; ky < obs_ensemble[0].Sq2D[kx].size(); ky++)
            {
                f << "," << avg_Sq2D_af[kx][ky];
            }
        }

        for (int kx = 0; kx < obs_ensemble[0].Sq2D.size(); kx++)
        {
            f << "\nSqSq2D,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            for (int ky = 0; ky < obs_ensemble[0].Sq2D[kx].size(); ky++)
            {
                f << "," << avg_SqSq2D[kx][ky];
            }
        }
        f.close();
    }
}

void ideal_gas::run_simulation(int N_config, int bin_num, std::string folder, std::string finfo)
{
    // run the simulation
    std::vector<observable> obs_ensemble(N_config);
    for (int i = 0; i < N_config; i++)
    {
        double percentage = (static_cast<double>(i + 1) / N_config) * 100;
        std::cout << "Generating configuration " << i + 1 << " of " << N_config << " (" << percentage << "%)\r" << std::flush;
        generate_gas();
        obs_ensemble[i] = measure_observable(bin_num);
    }
    save_gas_config_to_file(folder + "/config_" + finfo + ".csv");
    save_observable_to_file(folder + "/obs_" + finfo + ".csv", obs_ensemble);
}