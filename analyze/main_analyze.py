#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np
import os


def main():
    # test_plot()

    if 1:
        folder = "../data/data_local/data_pool"
        n, sigma, sqrtD, gxy = 100, 0.0, 1.0, 5.0
        parameters = [
            [100, 0.0, 1.0, 10.0],
            [150, 0.0, 1.0, 10.0],
        ]
        for parameter in parameters:
            n, sigma, sqrtD, gxy = parameter
            finfo = f"n{n:.0f}_sigma{sigma:.1f}_sqrtD{sqrtD:.1f}_gxy{gxy:.1f}"
            filename = folder + f"/config_{finfo}.csv"
            plot_gas_config(filename, finfo, show=True)
            # plot_gas_Sq_SqSq(folder, parameter, show=True)
            plot_gas_Iq_IqIq(folder, finfo, show=True)
    if 0:
        folder = "../data/data_local/data_pool"
        rand_num = 2
        for rnum in range(rand_num):
            finfo = f"random_run{rnum}"
            filename = folder + f"/config_{finfo}.csv"
            plot_gas_config(filename, finfo, show=True)
            # plot_gas_Sq_SqSq(folder, parameter, show=True)
            plot_gas_Iq_IqIq(folder, finfo, show=True)


if __name__ == "__main__":
    main()
