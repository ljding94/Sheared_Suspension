#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np
from ML_analyze import *
import sys
import random
import time


def main():

    print("analyzing data using ML model")
    folder = "../data/20241206_andes_rand"
    #folder = "../data/20241125_rand"
    rand_num = 6500
    rand_max = 6000
    parameters = []
    for i in range(rand_num):
        filename = f"{folder}/obs_random_run{i}.csv"
        if os.path.exists(filename):
            parameters.append([i])
        if len(parameters) >= rand_max:
            break
    print("parameters", parameters)
    print("total number of parameters", len(parameters))

    #plot_pddf_acf(folder, parameters, max_z=3.0, n_bin=30)

    calc_svd(folder, parameters)

    return 0
    random.shuffle(parameters)
    parameters_train = parameters[: int(0.7 * len(parameters))]
    parameters_test = parameters[int(0.7 * len(parameters)) :]

    #all_feature_mean, all_feature_std, all_gp_per_feature = GaussianProcess_optimization(folder, parameters_train)
    all_feature_names, all_feature_mean, all_feature_std, all_gp_per_feature = read_gp_and_feature_stats(folder)

    GaussianProcess_prediction(folder, parameters_test, all_feature_mean, all_feature_std, all_gp_per_feature)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
