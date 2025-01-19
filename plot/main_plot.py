#!/opt/homebrew/bin/python3
# from obs_plot import *
from Iq_gq_plot import *
from ML_plot import *


def main():

    plot_config_Iq()

    plot_Iq_gq_single()
    plot_Iq_gq_avg()

    plot_theory_gq()

    plot_SVD_data()
    plot_SVD_feature_data()
    plot_LML_contour()
    plot_GPR_prediction()


# todo: 1. increas hspace between subplot
# 2. remove legend for theory line
# 3. use Q instead of q
# 4. add kappa=10 in Sq plot
# 5. wor on TOC, use exact size as instructed


if __name__ == "__main__":
    main()
