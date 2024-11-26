import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
import os
from scipy.optimize import curve_fit
import pickle


def get_feature_Iq_gq_data(folder, parameters, random=True):

    all_Rmu = []  # system size related
    all_n = []  # system size related
    all_sigma = []
    all_sqrtD = []
    all_gxy = []

    all_Iq2D = []
    all_IqIq_af = []
    all_gq = []
    all_gq_phi = []
    all_gq_r = []

    qr = []
    qphi = []
    for i in range(len(parameters)):
        if random:
            filename = f"{folder}/obs_random_run{parameters[i][0]}.csv"
        else:
            N, sigma, theta, Sx, phi = parameters[i]
            finfo = f"N{N:.0f}_sigma{sigma:.1f}_theta{theta:.2f}_Sx{Sx:.2f}_phi{phi:.2f}"
            filename = f"{folder}/obs_{finfo}.csv"

        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue

        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        bnum_phi = len(data[0]) - 7
        print(f"bnum_phi: {bnum_phi}")
        bnum_r = int((len(data) - 2) / 3)
        print(f"bnum_r: {bnum_r}")

        Rmu, n, sigma, sqrtD, gxy = data[0, 1], data[0, 2], data[0, 3], data[0, 4], data[0, 5]

        qphi = data[1, 7:]
        qr = data[2 : 2 + bnum_r, 6]
        Iq2D = data[2 : 2 + bnum_r, 7:]
        Iq2D_af = data[2 + bnum_r : 2 + 2 * bnum_r, 7:]
        IqIq_af = data[2 + 2 * bnum_r : 2 + 3 * bnum_r, 7:]

        gq = IqIq_af / (Iq2D * Iq2D)
        QPHI, QR = np.meshgrid(qphi, qr)
        gq_phi = np.mean(gq*QR, axis=0)/np.mean(QR, axis=0)
        gq_r = np.mean(gq, axis=1)


        all_Iq2D.append(Iq2D)
        all_IqIq_af.append(IqIq_af)
        all_gq.append(gq)
        all_gq_phi.append(gq_phi)
        all_gq_r.append(gq_r)

        all_Rmu.append(Rmu)
        all_n.append(n)
        all_sigma.append(sigma)
        all_sqrtD.append(sqrtD)
        all_gxy.append(gxy)

    all_feature = np.array([all_Rmu, all_n, all_sigma, all_sqrtD, all_gxy]).T
    all_feature_name = ["Rmu", "n", "sigma", "sqrtD", "gxy"]
    all_feature_tex = [r"$R_\mu/R_0$", r"$n$", r"$R_\sigma$", r"$\sqrt{D}/R_0$", r"$\gamma_{xy}/R_0$"]

    qphi = np.array(qphi)
    qr = np.array(qr)
    all_Iq2D = np.array(all_Iq2D)
    all_IqIq_af = np.array(all_IqIq_af)
    all_gq = np.array(all_gq)
    all_gq_phi = np.array(all_gq_phi)
    all_gq_r = np.array(all_gq_r)

    return all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi

def calc_svd(folder, parameters):
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, parameters, random=True)

    #F = [np.concatenate((all_gq_r[i], all_gq_phi[i])) for i in range(len(all_gq_r))]
    F = all_gq_phi

    print("all_feature shape:", np.array(all_feature).shape)
    svd = np.linalg.svd(F)
    print(svd.S)
    print("np.array(svd.U).shape", np.array(svd.U).shape)
    print("np.array(svd.S).shape", np.array(svd.S).shape)
    print("np.array(svd.Vh).shape", np.array(svd.Vh).shape)
    # print(np.linalg.svd(all_Delta_Sq))
    plt.figure(figsize=(6, 6))
    # Subplot for svd.S
    ax00 = plt.subplot(2, 2, 1)
    ax00.plot(range(len(svd.S)), svd.S, "o--", markerfacecolor="none", label="svd.S")
    ax00.set_title("Singular Values (svd.S)")

    # Subplot for svd.U
    ax01 = plt.subplot(2, 2, 2)
    print("np.minimum(svd.Vh[0]), np.maximum(svd.Vh[0])", svd.Vh[0].min(), svd.Vh[0].max())
    #ax01.contourf(QPHI, QR, svd.Vh[0].reshape(np.shape(QPHI)), cmap="rainbow")  # , levels=np.linspace(-0.3, 0.2, 10), cmap="rainbow")
    ax01.plot(qphi, svd.Vh[0])
    ax01.set_title("Left Singular Vectors (svd.Vh[0])")

    ax11 = plt.subplot(2, 2, 3)
    print("np.minimum(svd.Vh[1]), np.maximum(svd.Vh[1])", svd.Vh[1].min(), svd.Vh[1].max())
    #ax11.contourf(QPHI, QR, svd.Vh[1].reshape(np.shape(QPHI)),cmap="rainbow")  # , levels=np.linspace(-0.3, 0.2, 10), cmap="rainbow")
    ax11.plot(qphi, svd.Vh[1])
    ax11.set_title("Left Singular Vectors (svd.Vh[1])")

    ax12 = plt.subplot(2, 2, 4)
    print("np.minimum(svd.Vh[2]), np.maximum(svd.Vh[2])", svd.Vh[2].min(), svd.Vh[2].max())
    ax12.plot()
    #ax12.contourf(QPHI, QR, svd.Vh[2].reshape(np.shape(QPHI)),cmap="rainbow") #, levels=np.linspace(-0.3, 0.2, 10), cmap="rainbow")
    ax12.plot(qphi, svd.Vh[2])
    ax12.set_title("Left Singular Vectors (svd.Vh[2])")

    plt.tight_layout()
    plt.savefig(f"{folder}/svd.png", dpi=300)
    plt.show()
    plt.close()

    SqV = np.inner(F, np.transpose(svd.Vh))
    plt.figure()
    fig = plt.figure(figsize=(2 * len(all_feature_name), 8))
    axs = [fig.add_subplot(2, len(all_feature_name) // 2 + 1, i + 1, projection="3d") for i in range(len(all_feature_name))]
    for i in range(len(all_feature_name)):
        scatter = axs[i].scatter(SqV[:, 0], SqV[:, 1], SqV[:, 2], c=all_feature[:, i], cmap="rainbow", s=1)
        axs[i].set_xlabel("V[0]")
        axs[i].set_ylabel("V[1]")
        axs[i].set_zlabel("V[2]")
        axs[i].set_title(all_feature_tex[i])
        axs[i].set_box_aspect([1, 1, 1])  # Set the aspect ratio of the plot
        # Set the same range for each axis
        max_range = np.array([SqV[:, 0].max() - SqV[:, 0].min(), SqV[:, 1].max() - SqV[:, 1].min(), SqV[:, 2].max() - SqV[:, 2].min()]).max() / 2.0
        mid_x = (SqV[:, 0].max() + SqV[:, 0].min()) * 0.5
        mid_y = (SqV[:, 1].max() + SqV[:, 1].min()) * 0.5
        mid_z = (SqV[:, 2].max() + SqV[:, 2].min()) * 0.5
        axs[i].set_xlim(mid_x - max_range, mid_x + max_range)
        axs[i].set_ylim(mid_y - max_range, mid_y + max_range)
        axs[i].set_zlim(mid_z - max_range, mid_z + max_range)
        cbar = fig.colorbar(scatter, ax=axs[i], fraction=0.02)
        cbar.set_label(all_feature_tex[i])
        axs[i].view_init(elev=10.0, azim=-30)

    plt.tight_layout()
    plt.savefig(f"{folder}/svd_projection_scatter_plot.png", dpi=300)
    plt.show()
    plt.close()

     # svd data
    data = np.column_stack((qphi, svd.S, svd.Vh[0], svd.Vh[1], svd.Vh[2]))
    column_names = ['q', 'svd.S', 'svd.Vh[0]', 'svd.Vh[1]', 'svd.Vh[2]']
    np.savetxt(f"{folder}/data_svd.txt", data, delimiter=',', header=','.join(column_names), comments='')

    data = np.column_stack((all_feature, SqV[:, 0], SqV[:, 1], SqV[:, 2]))
    column_names = all_feature_name + ["sqv[0]", "sqv[1]", "sqv[2]"]
    np.savetxt(f"{folder}/data_svd_projection.txt", data, delimiter=",", header=",".join(column_names), comments="")


def calc_Sq_pair_distance_distribution(all_Delta_Sq, max_z, bin_num):
    all_z = np.linspace(0, max_z, bin_num)
    all_Delta_Sq_dis = np.zeros(bin_num)
    all_Delta_Sq_dis[0] = 1 / len(all_Delta_Sq)  # for self distance

    for i in range(len(all_Delta_Sq) - 1):
        for j in range(i + 1, len(all_Delta_Sq)):
            Delta_Sq_dis = np.sqrt(np.sum(np.square(all_Delta_Sq[i] - all_Delta_Sq[j])))
            bin_index = int(Delta_Sq_dis / (max_z / bin_num))
            if bin_index >= bin_num:
                raise ValueError(f"bin_index >= bin_num, z={bin_index/bin_num*max_z}")
            all_Delta_Sq_dis[bin_index] += 2.0 / (len(all_Delta_Sq)) ** 2 / (max_z / bin_num)  # 2.0 for i,j and j,i symmetry, normalize to 1

    return all_Delta_Sq_dis, all_z


def calc_Sq_autocorrelation(mu, all_Delta_Sq, max_z, bin_num):
    # measure the autocorrelation of mu(Delta_Sq)
    all_z = np.linspace(0, max_z, bin_num)
    avg_mu = np.mean(mu)
    avg_mu2 = np.mean(np.square(mu))
    print("np.shape(mu)", np.shape(mu))
    print("np.shape(all_Delta_Sq)", np.shape(all_Delta_Sq))

    print("avg_mu:", avg_mu)
    print("avg_mu2:", avg_mu2)
    print("avg_mu2-avg_mu**2:", avg_mu2 - avg_mu**2)

    avg_mumuz = np.zeros(bin_num)
    avg_mu2z = np.zeros(bin_num)
    avg_muz = np.zeros(bin_num)

    # avg_mumuz[0] = avg_mu2  # for self distance
    # for i in range(len(all_Delta_Sq)-1):
    #    for j in range(i+1, len(all_Delta_Sq)):
    bin_count = np.zeros(bin_num)
    for i in range(len(mu)):
        for j in range(len(mu)):
            Delta_Sq_dis = np.sqrt(np.sum(np.square(all_Delta_Sq[i] - all_Delta_Sq[j])))
            bin_index = int(Delta_Sq_dis / (max_z / bin_num))
            if bin_index >= bin_num:
                raise ValueError(f"bin_index >= bin_num, z={bin_index/bin_num*max_z}")
            avg_muz[bin_index] += mu[i]
            avg_mu2z[bin_index] += mu[i] * mu[i]
            avg_mumuz[bin_index] += mu[i] * mu[j]
            bin_count[bin_index] += 1
    for i in range(bin_num):
        avg_muz[i] /= bin_count[i]
        avg_mu2z[i] /= bin_count[i]
        avg_mumuz[i] /= bin_count[i]

    ac_mu = np.ones(bin_num)
    for i in range(0, bin_num):
        if avg_mu2 - avg_mu**2 == 0:
            ac_mu[i] = 1
        else:
            ac_mu[i] = (avg_mumuz[i] - avg_muz[i] ** 2) / (avg_mu2z[i] - avg_muz[i] ** 2)

    if ac_mu[0] != 1:
        print("ac_mu[0]!=1: ")
        print("avg_mumuz[0]-avg_muz[0]**2,", avg_mumuz[0] - avg_muz[0] ** 2)
        print("(avg_mu2z[0]-avg_muz[0]**2)", (avg_mu2z[0] - avg_muz[0] ** 2))
        print(bin_count)
    ac_mu[0] = 1
    print("ac_mu", ac_mu)
    return ac_mu, all_z


def plot_pddf_acf(folder, parameters, max_z=2, n_bin=100):

    all_feature, all_feature_name, all_SqSq2D_flatten, qD = get_all_feature_Iq_gq_data(folder, parameters)

    p_z, z = calc_Sq_pair_distance_distribution(all_SqSq2D_flatten, max_z, n_bin)

    plt.figure(figsize=(8, 6))
    plt.plot(z, p_z / np.max(p_z), label="p_z/max(p_z)")

    acf_data = []
    for i in range(len(all_feature_name)):
        # pass
        acf_mu, z = calc_Sq_autocorrelation(all_feature[:, i], all_SqSq2D_flatten, max_z, n_bin)
        plt.plot(z, acf_mu, label=f"acf_{all_feature_name[i]}")
        acf_data.append(acf_mu)

    plt.xlabel("z")
    plt.ylabel("Value")
    plt.title("Pair Distance Distribution and Autocorrelation")
    plt.legend()
    plt.savefig(f"{folder}/pddf_acf.png", dpi=300)
    plt.close()

    # save these data to file for futher easy plotting

    data = np.column_stack((z, p_z, *acf_data))
    column_names = ["z", "p_z", *["acf_" + feature_name for feature_name in all_feature_name]]
    np.savetxt(f"{folder}/data_pddf_acf.txt", data, delimiter=",", header=",".join(column_names), comments="")


def GaussianProcess_optimization(folder, parameters_train):
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, parameters_train, random=True)

    grid_size = 30

    theta_per_feature = {
        "sigma": (np.logspace(-1, 0, grid_size), np.logspace(-3, -1, grid_size)),
        "sqrtD": (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
        "gxy": (np.logspace(-1, 1, grid_size), np.logspace(-4, -2, grid_size)),
    }

    # feature normalization
    all_feature_mean = np.mean(all_feature, axis=0)
    all_feature_std = np.std(all_feature, axis=0)
    all_feature = (all_feature - all_feature_mean) / all_feature_std
    all_gp_per_feature = {}
    plt.figure()
    fig, axs = plt.subplots(1, len(all_feature_name), figsize=(6 * len(all_feature_name), 6))
    for feature_name, (theta0, theta1) in theta_per_feature.items():
        if feature_name not in all_feature_name:
            continue
        print("training: ", feature_name)
        feature_index = all_feature_name.index(feature_name)

        F_learn = all_gq_phi

        # witout theta optimization
        kernel = RBF(1) + WhiteKernel(1)
        gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, optimizer=None).fit(F_learn, all_feature[:, feature_index])
        # print(" all_feature[:, feature_index]", all_feature[:, feature_index])

        print("GPML kernel: %s" % gp.kernel_)
        gp_theta = np.exp(gp.kernel_.theta)
        # kernel_params_array = np.array(list(kernel_params.values()))
        print("Kernel parameters:", gp_theta)
        print("Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))

        # calc Log likelihood
        ax = axs[all_feature_name.index(feature_name)]
        Theta0, Theta1 = np.meshgrid(theta0, theta1)
        LML = [[0 for j in range(Theta0.shape[1])] for i in range(Theta0.shape[0])]
        for i in range(Theta0.shape[0]):
            for j in range(Theta0.shape[1]):
                LML[i][j] = gp.log_marginal_likelihood(np.log([Theta0[i, j], Theta1[i, j]]))
                print(f"Calculating LML: i={i}/{Theta0.shape[0]}, j={j}/{Theta0.shape[1]}, LML={LML[i][j]}", end="\r")
        # reason for np.log here is the theta is log-transformed hyperparameters (https://github.com/scikit-learn/scikit-learn/blob/5491dc695/sklearn/gaussian_process/kernels.py#L1531) line (289)

        ax.contour(Theta0, Theta1, LML, levels=200)
        # find optimized theta0, theta1, using the above contour as guidanve
        kernel = RBF(theta0[grid_size // 2], (theta0[0], theta0[-1])) + WhiteKernel(theta1[grid_size // 2], (theta1[0], theta1[-1]))
        gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, n_restarts_optimizer=10).fit(F_learn, all_feature[:, feature_index])
        all_gp_per_feature[feature_name] = gp

        print("GPML kernel: %s" % gp.kernel_)
        gp_theta = np.exp(gp.kernel_.theta)
        # kernel_params_array = np.array(list(kernel_params.values()))
        print("Kernel parameters:", gp_theta)
        print("Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))

        ax.plot([gp_theta[0]], [gp_theta[1]], "x", color="red", markersize=10, markeredgewidth=2, label=r"l=%.2e, $\sigma$=%.2e" % (gp_theta[0], gp_theta[1]))

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"theta0: l")
        ax.set_ylabel(r"theta1: $\sigma$")
        feature_name_legend = feature_name
        ax.set_title(f"Log Marginal Likelihood for {feature_name_legend}")
        ax.legend()

        data = np.column_stack(([gp_theta[0]] * len(theta0), [gp_theta[1]] * len(theta1), theta0, theta1, np.array(LML).T))
        column_names = ["gp_theta0", "gp_theta1", "theta0", "theta1", "LML"]
        np.savetxt(f"{folder}/data_{feature_name}_LML.txt", data, delimiter=",", header=",".join(column_names), comments="")
        with open(f"{folder}/gp_{feature_name}.pkl", "wb") as f:
            pickle.dump(gp, f)

    # Save average and standard deviation per feature
    avg_std_data = np.column_stack((all_feature_name, all_feature_mean, all_feature_std))
    column_names = ["Feature", "Mean", "Std"]
    np.savetxt(f"{folder}/data_feature_avg_std.txt", avg_std_data, delimiter=",", header=",".join(column_names), comments="", fmt="%s")

    plt.tight_layout()
    plt.savefig(f"{folder}/LML_subplots.png", dpi=300)
    # plt.show()
    plt.close()

    # return trained GPR
    return all_feature_mean, all_feature_std, all_gp_per_feature


def read_gp_and_feature_stats(folder):
    all_feature_name = ["Rmu", "n", "sigma", "sqrtD", "gxy"]
    all_feature_tex = [r"$R_\mu/R_0$", r"$n$", r"$R_\sigma$", r"$\sqrt{D}/R_0$", r"$\gamma_{xy}/R_0$"]
    all_feature_mean = np.genfromtxt(f"{folder}/data_feature_avg_std.txt", delimiter=",", skip_header=1, usecols=1)
    all_feature_std = np.genfromtxt(f"{folder}/data_feature_avg_std.txt", delimiter=",", skip_header=1, usecols=2)
    all_gp_per_feature = {}
    for feature_name in all_feature_name:
        if os.path.exists(f"{folder}/gp_{feature_name}.pkl"):
            with open(f"{folder}/gp_{feature_name}.pkl", "rb") as f:
                all_gp_per_feature[feature_name] = pickle.load(f)
    return all_feature_name, all_feature_mean, all_feature_std, all_gp_per_feature


def GaussianProcess_prediction(folder, parameters_test, all_feature_mean, all_feature_std, all_gp_per_feature):
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, parameters_test, random=True)

    plt.figure()
    Fstar = all_gq_phi
    fig, axs = plt.subplots(1, len(all_feature_name), figsize=(6 * len(all_feature_name), 6))
    for feature_name, gp in all_gp_per_feature.items():
        feature_index = all_feature_name.index(feature_name)
        Y = all_feature[:, feature_index]

        print("GPML kernel: %s" % gp.kernel_)
        gp_theta = np.exp(gp.kernel_.theta)  # gp.kernel_.theta return log transformed theta
        # kernel_params_array = np.array(list(kernel_params.values()))
        print("Kernel parameters:", gp_theta)
        print("Log-marginal-likelihood: %.3f" % gp.log_marginal_likelihood(gp.kernel_.theta))

        Y_predict, Y_predict_err = gp.predict(Fstar, return_std=True)
        # print("np.shape(test_data[:, 0])", np.shape(test_data[:, 0]))
        print("np.shape(Fstar)", np.shape(Fstar))
        print("np.shape(Y_predict)", np.shape(Y_predict))

        Y_predict = Y_predict * all_feature_std[feature_index] + all_feature_mean[feature_index]
        Y_predict_err = Y_predict_err * all_feature_std[feature_index]

        axs[feature_index].errorbar(Y, Y_predict, yerr=Y_predict_err, marker="o", markerfacecolor="none", markersize=3, linestyle="none")
        axs[feature_index].plot(Y, Y, "--")
        min_val = min(np.min(Y), np.min(Y_predict - Y_predict_err))
        max_val = max(np.max(Y), np.max(Y_predict + Y_predict_err))
        axs[feature_index].set_xlim(min_val, max_val)
        axs[feature_index].set_ylim(min_val, max_val)
        axs[feature_index].set_xlabel(f"{feature_name}")
        axs[feature_index].set_ylabel(f"{feature_name} Prediction")

        # save data to file
        data = np.column_stack((Y, Y_predict, Y_predict_err))
        column_names = [feature_name, "ML predicted", "ML predicted uncertainty"]
        np.savetxt(f"{folder}/data_{feature_name}_prediction.txt", data, delimiter=",", header=",".join(column_names), comments="")

    plt.savefig(f"{folder}/prediction.png", dpi=300)
    plt.close()


def ax_fit(x, a):
    return a * x


def fit_Rg2(q, Sq):
    popt, pcov = curve_fit(ax_fit, q**2 / 3, (1 - Sq))
    perr = np.sqrt(np.diag(pcov))
    return popt[0], perr[0]


def calc_Sq_fitted_Rg2(folder, parameters_test, all_feature_name):
    segment_type, all_feature, all_feature_name, all_Sq, all_Sq_err, q = read_Sq_data(folder, parameters_test)

    MC_Rg2 = all_feature[:, all_feature_name.index("Rg2")]
    # qfns = [10,20,30,40]
    qfns = [50, 55, 60, 65, 70]
    Rg2s = []
    Rg2_errs = []
    plt.figure()
    for qfn in qfns:
        Rg2s.append([])
        Rg2_errs.append([])
        for i in range(len(all_Sq)):
            Rg2, Rg2_err = fit_Rg2(q[:qfn], all_Sq[i][:qfn])
            Rg2s[-1].append(Rg2)
            Rg2_errs[-1].append(Rg2_err)

        plt.scatter(MC_Rg2, Rg2s[-1], alpha=0.5, label=f"qf={q[qfn-1]}")
    plt.plot(MC_Rg2, MC_Rg2, "k--")
    plt.xlabel("MC Rg2")
    plt.ylabel("Fitted Rg2")
    plt.legend()
    plt.savefig(f"{folder}/{segment_type}_Rg2_fit.png", dpi=300)
    plt.close()

    data = np.column_stack(([MC_Rg2] + Rg2s))
    column_names = ["MC Rg2", "fitted Rg2"]
    np.savetxt(f"{folder}/data_{segment_type}_fitted_Rg2.txt", data, delimiter=",", header=",".join(column_names), comments="")
