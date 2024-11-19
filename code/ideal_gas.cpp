#include "ideal_gas.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>

// initialization
ideal_gas::ideal_gas(int n_, double R0_, double sqrtD0_, double sigma_, double gxx_, double gxy_, double gyx_, double gyy_, bool random_param)
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
    R0 = R0_;
    sqrtD0 = sqrtD0_;

    // generate the gas system
    if (random_param)
    {
        std::cout << "\nrandom_param" << std::endl;
        // sigma = 1.0 * rand_uni(gen);
        // theta = 0.1 * M_PI * rand_uni(gen);
        // Sx = 0.9 + 0.2 * rand_uni(gen);
        // phi = 0.1 * M_PI * rand_uni(gen);
        //  to make area invariant
        // gxx*gyy-gxy*gyx=1
        sigma = 0.0 * rand_uni(gen);
        gxx = 0.9 + 0.2 * rand_uni(gen);
        gxy = -0.1 + 0.2 * rand_uni(gen);
        gyx = -0.1 + 0.2 * rand_uni(gen);
        gyy = (1 + gxy * gyx) / gxx;
    }
    else
    {
        gxx = gxx_;
        gxy = gxy_;
        gyx = gyx_;
        gyy = gyy_;
        sigma = sigma_;
        // theta = theta_;
        // Sx = Sx_;
        // phi = phi_;
    }
    // Sy = 1.0 / Sx;
    // std::cout << "theta: " << theta << ", Sx: " << Sx << ", Sy: " << Sy << ", phi: " << phi << std::endl;
    std::cout << "\ngxx: " << gxx << ", gxy: " << gxy << ", gyx: " << gyx << ", gyy: " << gyy << std::endl;

    // calc the affine matrix
    /*
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);
    */
    // gxx = Sx * cos_theta * cos_phi + Sy * sin_theta * sin_phi;
    // gxy = Sx * sin_theta * cos_phi - Sy * cos_theta * sin_phi;
    // gyx = Sx * cos_theta * sin_phi - Sy * sin_theta * cos_phi;
    // gyy = Sx * sin_theta * sin_phi + Sy * cos_theta * cos_phi;
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
        all_beads[i].R = R0 * std::exp(sigma * rand_norm(gen)); // on average radius of the bead PDI=e^{3*sigma^2}: I guess
        all_beads[i].V = 4.0 / 3.0 * M_PI * all_beads[i].R * all_beads[i].R * all_beads[i].R;
        all_beads[i].r = {(rand_uni(gen) - 0.5) * L, (rand_uni(gen) - 0.5) * L}; // bead position in 2d
        //  let's try circular shape
        // r_r = L * std::sqrt(rand_uni(gen));
        // r_theta = 2 * M_PI * rand_uni(gen);
        // all_beads[i].r = {r_r * std::cos(r_theta), r_r * std::sin(r_theta)};
    }
    return 1;
}

observable ideal_gas::measure_observable(int bnum_r, int bnum_phi)
{
    // measure the observable
    observable obs;

    // measure the structure factor

    // 1. set up q vectors
    double qri = 500;
    double qrf = 5000;
    obs.qr.resize(bnum_r);
    std::vector<double> qr(bnum_r, 0);
    for (int k = 0; k < bnum_r; k++)
    {
        obs.qr[k] = qri * std::pow(qrf / qri, 1.0 * k / (bnum_r - 1)); // uniform in log scale
        // obs.qr[k] = qri + (qrf - qri) * k / (bnum_r - 1); // uniform in linear scale
    }
    obs.qphi.resize(bnum_phi);
    for (int k = 0; k < bnum_phi; k++)
    {
        obs.qphi[k] = M_PI * k / (bnum_phi - 1);
    }

    // 2. calculate form factor for each bead
    for (int i = 0; i < all_beads.size(); i++)
    {
        all_beads[i].fq = calc_structure_sphere_form_factor(obs.qr, all_beads[i].R);
    }
    obs.Iq2D = calc_structure_factor_2d(obs.qr, obs.qphi);
    obs.Iq2D_af = calc_structure_factor_2d_af(obs.qr, obs.qphi);

    std::vector<std::vector<double>> IqIq(bnum_r, std::vector<double>(bnum_phi, 0.0));
    std::vector<std::vector<double>> Iq_afIq_af(bnum_r, std::vector<double>(bnum_phi, 0.0));
    std::vector<std::vector<double>> IqIq_af(bnum_r, std::vector<double>(bnum_phi, 0.0));

    IqIq = element_wise_dot_2d(obs.Iq2D, obs.Iq2D);
    Iq_afIq_af = element_wise_dot_2d(obs.Iq2D_af, obs.Iq2D_af);
    IqIq_af = element_wise_dot_2d(obs.Iq2D, obs.Iq2D_af);

    // ref: https://en.wikipedia.org/wiki/Dynamic_light_scattering

    /*
    for (int kr = 0; kr < bnum_r; kr++)
    {
        for (int kphi = 0; kphi < bnum_phi; kphi++)
        {
            // IqIq_af[kr][kphi] /= 1.0*IqIq[kr]d[kphi];
            // IqIq_af[kr][kphi] /= 0.5*(IqIq[kr][kphi] + Iq_afIq_af[kr][kphi]);
            IqIq_af[kr][kphi] /= std::pow(obs.Iq2D[kr][kphi] + obs.Iq2D_af[kr][kphi], 2) / 4;
        }
    }
    */

    obs.IqIq_af = IqIq_af;

    return obs;
}

std::vector<double> ideal_gas::calc_structure_sphere_form_factor(std::vector<double> qr, double R)
{
    int bnum_r = qr.size();
    std::vector<double> fq(bnum_r, 0);
    for (int k = 0; k < bnum_r; k++)
    {
        fq[k] = 3.0 / std::pow((qr[k] * R), 3) * (std::sin(qr[k] * R) - qr[k] * R * std::cos(qr[k] * R));
    }
    return fq;
}
std::vector<std::vector<double>> ideal_gas::calc_structure_factor_2d(std::vector<double> qr, std::vector<double> qphi)
{
    // measrue 2d structure factor
    int bnum_r = qr.size();
    int bnum_phi = qphi.size();

    std::vector<double> r{0, 0};
    double qx, qy, Aq_Re_buff, Aq_Im_buff;
    std::vector<std::vector<double>> Iq(bnum_r, std::vector<double>(bnum_phi, 0.0));
    std::vector<std::vector<double>> Aq_Re(bnum_r, std::vector<double>(bnum_phi, 0.0));
    std::vector<std::vector<double>> Aq_Im(bnum_r, std::vector<double>(bnum_phi, 0.0));

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

    // calculate Sq

    double sum_V2 = 0; // sum of V^2 of each scatter bead, for final normalization
    for (int i = 0; i < scatter_beads.size(); i++)
    {
        sum_V2 += scatter_beads[i].V * scatter_beads[i].V;
    }

    for (int kr = 0; kr < bnum_r; kr++)
    {
        for (int kphi = 0; kphi < bnum_phi; kphi++)
        {

            qx = qr[kr] * std::cos(qphi[kphi]);
            qy = qr[kr] * std::sin(qphi[kphi]);
            for (int i = 0; i < scatter_beads.size(); i++)
            {
                r[0] = scatter_beads[i].r[0];
                r[1] = scatter_beads[i].r[1];
                Aq_Re_buff = scatter_beads[i].V * scatter_beads[i].fq[kr] * std::cos(qx * r[0] + qy * r[1]);
                Aq_Re[kr][kphi] += Aq_Re_buff;
                Aq_Im_buff = scatter_beads[i].V * scatter_beads[i].fq[kr] * std::sin(qx * r[0] + qy * r[1]);
                Aq_Im[kr][kphi] += Aq_Im_buff;
            }
            Iq[kr][kphi] = Aq_Re[kr][kphi] * Aq_Re[kr][kphi] + Aq_Im[kr][kphi] * Aq_Im[kr][kphi];
            Iq[kr][kphi] /= sum_V2;
        }
    }
    return Iq;
}

void ideal_gas::affine_transform()
{
    // affine transform the gas system
    double x, y;
    for (int i = 0; i < all_beads.size(); i++)
    {
        x = all_beads[i].r[0];
        y = all_beads[i].r[1];
        all_beads[i].r[0] = gxx * x + gxy * y;
        all_beads[i].r[1] = gyx * x + gyy * y;
    }
}

void ideal_gas::Brownian_transform()
{
    // Brownian transform the gas system
    double dx, dy;
    for (int i = 0; i < all_beads.size(); i++)
    {
        dx = sqrtD0 * std::sqrt(R0 / all_beads[i].R) * rand_norm(gen);
        dy = sqrtD0 * std::sqrt(R0 / all_beads[i].R) * rand_norm(gen);
        all_beads[i].r[0] += dx;
        all_beads[i].r[1] += dy;
    }
}

std::vector<std::vector<double>> ideal_gas::calc_structure_factor_2d_af(std::vector<double> qr, std::vector<double> qphi)
{
    // measrue 2d structure factor after affine transformation
    std::vector<std::vector<double>> Sq; // initialize with 1 due to self overlaping term (see S.H. Chen 1986 eq 18)
    std::vector<bead> all_beads_pre_af = all_beads;
    affine_transform();
    Brownian_transform();
    Sq = calc_structure_factor_2d(qr, qphi);
    all_beads = all_beads_pre_af;
    return Sq;
}

std::vector<std::vector<double>> ideal_gas::element_wise_dot_2d(std::vector<std::vector<double>> M1, std::vector<std::vector<double>> M2)
{
    // calculate the dot product of two 2d structure factors
    int bnum_r = M1.size();
    int bnum_phi = M1[0].size();
    std::vector<std::vector<double>> MM(bnum_r, std::vector<double>(bnum_phi, 0.0));
    for (int kr = 0; kr < bnum_r; kr++)
    {
        for (int kphi = 0; kphi < bnum_phi; kphi++)
        {
            MM[kr][kphi] = M1[kr][kphi] * M2[kr][kphi];
        }
    }
    return MM;
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
        Brownian_transform();
        f << "R,x,y,x_af,y_af\n";
        for (int i = 0; i < all_beads.size(); i++)
        {
            f << all_beads_pre_af[i].R << "," << all_beads_pre_af[i].r[0] << "," << all_beads_pre_af[i].r[1] << "," << all_beads[i].r[0] << "," << all_beads[i].r[1] << "\n";
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
        std::vector<std::vector<double>> avg_Iq2D(obs_ensemble[0].Iq2D.size(), std::vector<double>(obs_ensemble[0].Iq2D[0].size(), 0.0));
        std::vector<std::vector<double>> avg_Iq2D_af(obs_ensemble[0].Iq2D.size(), std::vector<double>(obs_ensemble[0].Iq2D[0].size(), 0.0));
        std::vector<std::vector<double>> avg_IqIq_af(obs_ensemble[0].Iq2D.size(), std::vector<double>(obs_ensemble[0].Iq2D[0].size(), 0.0));

        double M = obs_ensemble.size();
        double sqrt_M = std::sqrt(M);

        // get the average
        for (int i = 0; i < obs_ensemble.size(); i++)
        {
            for (int kr = 0; kr < obs_ensemble[i].Iq2D.size(); kr++)
            {
                for (int kphi = 0; kphi < obs_ensemble[i].Iq2D[kr].size(); kphi++)
                {
                    avg_Iq2D[kr][kphi] += obs_ensemble[i].Iq2D[kr][kphi];
                    avg_Iq2D_af[kr][kphi] += obs_ensemble[i].Iq2D_af[kr][kphi];
                    avg_IqIq_af[kr][kphi] += obs_ensemble[i].IqIq_af[kr][kphi];
                }
            }
        }

        for (int kr = 0; kr < obs_ensemble[0].Iq2D.size(); kr++)
        {
            for (int kphi = 0; kphi < obs_ensemble[0].Iq2D[kr].size(); kphi++)
            {
                avg_Iq2D[kr][kphi] /= M;
                avg_Iq2D_af[kr][kphi] /= M;
                avg_IqIq_af[kr][kphi] /= M;
            }
        }

        f << "label,n,sigma,theta,Sx,Sy,phi,gxx,gxy,gyz,gyy,qr,Iq2D/Iq2D_af/IqIq_af\n";
        // write parameters and stats to the file
        f << "mean," << n << "," << sigma << "," << theta << "," << Sx << "," << Sy << "," << phi
          << "," << gxx << "," << gxy << "," << gyx << "," << gyy;
        f << ",Na";
        for (int j = 0; j < obs_ensemble[0].qphi.size(); j++)
        {
            f << ",NA";
        }
        f << "\n qphi,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
        f << ",Na";
        for (int j = 0; j < obs_ensemble[0].qphi.size(); j++)
        {
            f << "," << obs_ensemble[0].qphi[j];
        }

        for (int kr = 0; kr < obs_ensemble[0].Iq2D.size(); kr++)
        {
            f << "\nIq2D,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            f << "," << obs_ensemble[0].qr[kr];
            for (int kphi = 0; kphi < obs_ensemble[0].Iq2D[kr].size(); kphi++)
            {
                f << "," << avg_Iq2D[kr][kphi];
            }
        }

        for (int kr = 0; kr < obs_ensemble[0].Iq2D_af.size(); kr++)
        {
            f << "\nIq2D_af,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            f << "," << obs_ensemble[0].qr[kr];
            for (int kphi = 0; kphi < obs_ensemble[0].Iq2D[kr].size(); kphi++)
            {
                f << "," << avg_Iq2D_af[kr][kphi];
            }
        }

        for (int kr = 0; kr < obs_ensemble[0].IqIq_af.size(); kr++)
        {
            f << "\nIqIq_af,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA";
            f << "," << obs_ensemble[0].qr[kr];
            for (int kphi = 0; kphi < obs_ensemble[0].Iq2D[kr].size(); kphi++)
            {
                f << "," << avg_IqIq_af[kr][kphi];
            }
        }
        f.close();
    }
}

void ideal_gas::run_simulation(int N_config, int bnum_r, int bnum_phi, std::string folder, std::string finfo)
{
    // run the simulation
    std::vector<observable> obs_ensemble(N_config);
    for (int i = 0; i < N_config; i++)
    {
        double percentage = (static_cast<double>(i + 1) / N_config) * 100;
        std::cout << "Generating configuration " << i + 1 << " of " << N_config << " (" << percentage << "%)\r" << std::flush;
        generate_gas();
        obs_ensemble[i] = measure_observable(bnum_r, bnum_phi);
    }
    save_gas_config_to_file(folder + "/config_" + finfo + ".csv");
    save_observable_to_file(folder + "/obs_" + finfo + ".csv", obs_ensemble);
}