#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np
import os


def main():
    folder = "../data/data_local/data_pool"
    N = 200
    parameters = [
        #[N, 0.0, 0.1, 0.9, 0.1],
        #[N, 0.0, 0.1, 0.9, 0.0],
        [N, 0.0, 0.1, 1.0, 0.0],
        [N, 0.0, 0.2, 1.0, 0.0],
        [N, 0.0, 0.0, 0.9, 0.0],
        [N, 0.0, 0.0, 0.8, 0.0],
        #[L, N, 1.0, 0.5, 0.5, 0.0],
        #[L, N, 0.0, 0.0, 1.0, 0.0],
        #[L, N, 0.0, 0.5, 1.0, 0.0],
        #[L, N, 1.0, 0.0, 1.0, 0.0],
        #[L, N, 1.0, 0.5, 1.0, 0.0],
        #[L, N, 0.5, 1.0, 0.0],
        #[L, N, 0.0, 2.0, 0.0],
        #[L, N, 0.0, 2.0, 0.5, 0.0],

        #[L, N, 0.0, 1.5, 1.0, 0.0],
        #[L, N, 0.0, 1.0, 1.5, 0.0],
        #[L, N, 0.0, 1.0, 1.0, 0.1],
        #[L, N, 0.0, 1.3, 1.0, 0.0],
    ]
    for parameter in parameters:
        N, sigma, theta, Sx, phi = parameter
        finfo = f"N{N:.0f}_sigma{sigma:.1f}_theta{theta:.1f}_Sx{Sx:.1f}_phi{phi:.1f}"
        filename = folder + f"/config_{finfo}.csv"
        plot_gas_config(filename, finfo, show=True)
        plot_gas_Sq_SqSq(folder, parameter, show=True)


if __name__ == "__main__":
    main()
