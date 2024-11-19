#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np
import os


def main():
    #test_plot()

    if 1:
        folder = "../data/data_local/data_pool"
        N = 20
        parameters = [
            [N, 0.0, 1.00, 0.00, 0.00, 1.00],
            [N, 0.0, 1.00, 0.10, 0.00, 1.00],
            [N, 0.0, 1.00, 0.20, 0.00, 1.00],
            #[N, 0.0, 1.00, 0.20, 0.00, 1.00],
            #[N, 0.0, 1.00, 0.30, 0.00, 1.00],
            #[N, 0.0, 1.00, 1.00, 0.00, 1.00],
            #[N, 0.1, 1.00, 0.10, 0.00, 1.00],
            #[N, 1.0, 1.00, 0.10, 0.00, 1.00],
        ]
        for parameter in parameters:
            N, sigma, gxx, gxy, gyx, gyy = parameter
            finfo = f"N{N:.0f}_sigma{sigma:.1f}_gxx{gxx:.2f}_gxy{gxy:.2f}_gyx{gyx:.2f}_gyy{gyy:.2f}"
            filename = folder + f"/config_{finfo}.csv"
            plot_gas_config(filename, finfo, show=True)
            # plot_gas_Sq_SqSq(folder, parameter, show=True)
            plot_gas_Iq_IqIq(folder, finfo, show=True)
    if 0:
        folder = "../data/data_local/data_pool"
        N = 200
        rand_num = 2
        for rnum in range(rand_num):
            finfo = f"N{N:.0f}_random_run{rnum}"
            filename = folder + f"/config_{finfo}.csv"
            plot_gas_config(filename, finfo, show=True)
            # plot_gas_Sq_SqSq(folder, parameter, show=True)
            plot_gas_Iq_IqIq(folder, finfo, show=True)


if __name__ == "__main__":
    main()
