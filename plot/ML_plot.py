from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import os
from Iq_gq_plot import get_feature_Iq_gq_data


def read_SVD_data(folder):
    data = np.loadtxt(f"{folder}/data_svd_projection.txt", skiprows=1, delimiter=",", unpack=True)
    print("data.shape", data.shape)
    print("data[:,-3:].shape", data[-3:, :].shape)
    sqv0, sqv1, sqv2 = data[-3:, :]

    all_feature = data[:-3, :]
    return all_feature, sqv0, sqv1, sqv2


def calc_sum_log_likelihood(x, k):
    mu1 = np.mean(x[:k])
    mu2 = np.mean(x[k:])
    sigma1 = np.std(x[:k])
    sigma2 = np.std(x[k:])
    q = len(x)
    sigma = np.sqrt(((k - 1) * sigma1**2 + (q - k - 1) * sigma2**2) / (q - 2))
    log_norm_x1 = stats.norm.logpdf(x[:k], mu1, sigma)
    log_norm_x2 = stats.norm.logpdf(x[k:], mu2, sigma)
    return np.sum(log_norm_x1) + np.sum(log_norm_x2)


def calc_maximum_profile_likelihood(svdS):
    k = range(1, 20)
    pll = np.zeros(len(k))
    for i in range(len(k)):
        pll[i] = calc_sum_log_likelihood(svdS, k[i])
    return k, pll


def plot_SVD_data(tex_lw=240.71031, ppi=72):
    #folder = "../data/20241125_rand"
    folder = "../data/20241206_andes_rand"
    data = np.loadtxt(f"{folder}/data_svd.txt", skiprows=1, delimiter=",", unpack=True)
    q, S, V0, V1, V2 = data
    #V0, V1, V2 = -V0, -V1, -V2

    plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.5))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((1, 2), (0, 0))
    ax01 = plt.subplot2grid((1, 2), (0, 1))
    # ax10 = plt.subplot2grid((2, 2), (1, 0), colspan=2)

    # Subplot for svd.S
    # ax00.plot(range(len(svd.S[:51])), svd.S[:51], "x", markersize=5, markerfacecolor='none', label=r"$\Sigma$")
    ax00.semilogx(range(1, len(S) + 1), S, "x", markersize=5, markerfacecolor="none", label=r"$\Sigma$")
    ax00.plot(range(1, len(S) + 1)[:3], S[:3], "ro", markersize=5, markerfacecolor="none")
    ax00.set_xlabel(r"$SVR$", fontsize=9)  # , labelpad=0)
    ax00.set_ylabel(r"$\Sigma$", fontsize=9, labelpad=0)
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax00.yaxis.set_major_locator(plt.MultipleLocator(200))
    ax00.yaxis.set_minor_locator(plt.MultipleLocator(100))

    ax01.plot(q, V0, "-", lw=1, label=r"$V0$")
    ax01.plot(q, V1, "-", lw=1, label=r"$V1$")
    ax01.plot(q, V2, "-", lw=1, label=r"$V2$")
    ax01.set_xlabel(r"$\theta$", fontsize=9)  # , labelpad=0)
    ax01.set_ylabel(r"$V$", fontsize=9, labelpad=0)
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax01.legend(loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.5, frameon=False, fontsize=9)
    ax01.yaxis.set_major_locator(plt.MultipleLocator(0.3))
    ax01.yaxis.set_minor_locator(plt.MultipleLocator(0.15))
    #ax01.set_ylim(-0.35, 0.4)
    phi_ticks = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
    phi_ticks_label = [r"$0$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"]
    ax01.set_xticks(phi_ticks)
    ax01.set_xticklabels(phi_ticks_label)
    """
    gfolder = "../data/20241130_prec"
    gfinfo = "n100_Rmu1.0_sigma0.0_sqrtD1.0_gxy10.0"
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(gfolder, [gfinfo])
    gphi = all_gq_phi[0]
    gphi0 = np.inner(gphi, V0) * V0
    gphi1 = np.inner(gphi, V1) * V1
    gphi2 = np.inner(gphi, V2) * V2

    ax10.plot(q, gphi, "-", color="red", lw=1, label=r"$\phi$")
    ax10.plot(q, gphi0+gphi1+gphi2, linestyle=(2, (2, 2)), color="blue", lw=1, label=r"$\phi$")
    ax10.plot(q, gphi0, "-", lw=1, label=r"$g0$")
    ax10.plot(q, gphi1, "-", lw=1, label=r"$g1$")
    ax10.plot(q, gphi2, "-", lw=1, label=r"$g2$")
    """

    ax00.text(0.8, 0.15, r"$(a)$", fontsize=9, transform=ax00.transAxes)
    ax01.text(0.8, 0.15, r"$(b)$", fontsize=9, transform=ax01.transAxes)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/SVD.pdf", format="pdf", dpi=300)
    plt.savefig("figures/SVD.png", dpi=300)

    plt.show()


def plot_SVD_feature_data(tex_lw=240.71031, ppi=72):
    #folder = "../data/20241125_rand"
    folder = "../data/20241206_andes_rand"
    all_feature, sqv0, sqv1, sqv2 = read_SVD_data(folder)
    #sqv0, sqv1, sqv2 = -sqv0, -sqv1, -sqv2
    all_feature_name = ["Rmu", "n", "sigma", "sqrtD", "gxy"]
    # all_feature_tex = [r"$R_\mu/R_0$", r"$n$", r"$R_\sigma$", r"$\sqrt{D}/R_0$", r"$\gamma/R_0$"]
    feature_to_tex = {"n": r"$nL^2$", "sigma": r"$R_s$", "sqrtD": r"$D_2$", "gxy": r"$\gamma L$"}

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    # n sigma, sqrtD, gxy

    ax1 = fig.add_subplot(221, projection="3d")
    ax2 = fig.add_subplot(222, projection="3d")
    ax3 = fig.add_subplot(223, projection="3d")
    ax4 = fig.add_subplot(224, projection="3d")


    '''
    ax1 = fig.add_subplot(141, projection="3d")
    ax2 = fig.add_subplot(142, projection="3d")
    ax3 = fig.add_subplot(143, projection="3d")
    ax4 = fig.add_subplot(144, projection="3d")
    '''
    axs = [ax1, ax2, ax3, ax4]

    cbar_lim=[(100,200),(0,0.3),(0.5,3),(5,30)]

    i = -1
    for feature_name, feature_tex in feature_to_tex.items():
        i += 1
        # cbar_major_locator = [0.4, 0.03, 5, 5, 40, 0.2, 0.1, 40]
        if feature_name not in all_feature_name:
            continue
        print("plotting: ", feature_name)

        feature_index = all_feature_name.index(feature_name)
        mu = all_feature[feature_index, :]

        ax = axs[i]
        ax.set_box_aspect((1, 1, 1.15))
        scatter = ax.scatter(sqv0, sqv1, sqv2, s=0.5, c=mu, cmap="rainbow", rasterized=True)

        ax.set_xlabel(r"$FV0$", fontsize=9, labelpad=-10, rotation=0)
        ax.set_ylabel(r"$FV1$", fontsize=9, labelpad=-11, rotation=0)
        ax.set_zlabel(r"$FV2$", fontsize=9, labelpad=-8, rotation=0, zorder=10)
        #ax.zaxis.set_label_coords(0.5, -0.5, 1.1)

        # ax.tick_params(labelsize=7, pad=0)
        ax.tick_params("x", labelsize=7, pad=-5)
        ax.tick_params("y", labelsize=7, pad=-5)
        ax.tick_params("z", labelsize=7, pad=-2)

        # ax.set_title(features_tex[i])
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.025, pad=-0.1, location="left")  # , location="top", orientation='horizontal')
        #cbar = fig.colorbar(scatter, ax=ax, fraction=0.025, pad=-0.1, location="right")  # , location="top", orientation='horizontal')
        cbar.ax.tick_params(direction="in", labelsize=7, pad=2)
        cbar.ax.set_ylim(cbar_lim[i])

        # cbar.ax.yaxis.set_major_locator(plt.MultipleLocator(cbar_major_locator[i]))
        # cbar.ax.yaxis.set_minor_locator(plt.MultipleLocator(cbar_major_locator[i]*0.5))

        # cbar.ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        # cbar.ax.set_axes_locator(plt.MultipleLocator(cbar_major_locator[i]))
        # cbar.ax.xaxis.set_minor_locator(plt.MultipleLocator(cbar_major_locator[i]*0.5))
        cbar.ax.set_title(feature_tex, fontsize=9)

        # cbar = .colorbar(axs.collections[0])
        # cbar.set_label(mu, fontsize=9)

        # ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        # ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
        # ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
        # ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        # ax.zaxis.set_major_locator(plt.MultipleLocator(0.4))
        # ax.zaxis.set_minor_locator(plt.MultipleLocator(0.2))
        #ax.view_init(elev=20.0, azim=124)
        ax.view_init(elev=17.0, azim=-55)
        ax.grid(True, which="minor")

    # plt.tight_layout( h_pad=0.2, w_pad=1.7) #, h_pad=-3, w_pad=2)

    ax1.text2D(0.8, 0.85, r"$(a)$", fontsize=9, transform=ax1.transAxes)
    ax2.text2D(0.8, 0.85, r"$(b)$", fontsize=9, transform=ax2.transAxes)
    ax3.text2D(0.8, 0.85, r"$(c)$", fontsize=9, transform=ax3.transAxes)
    ax4.text2D(0.8, 0.85, r"$(d)$", fontsize=9, transform=ax4.transAxes)

    #plt.tight_layout(h_pad=1, w_pad=3.2, pad=1)
    plt.subplots_adjust(left=0.075, bottom=0.05, right=0.955, top=0.95, wspace=0.2, hspace=0.2)
    plt.savefig("figures/SVD_feature.pdf", format="pdf", dpi=300, transparent=True)
    plt.savefig("figures/SVD_feature.png", dpi=300)
    plt.show()
    plt.close()


def get_LML_date(folder, feature):
    data = np.loadtxt(f"{folder}/data_{feature}_LML.txt", skiprows=1, delimiter=",", unpack=True)
    gp_theta0, gp_theta1, theta0, theta1, LML = data[0], data[1], data[2], data[3], data[4:]
    return gp_theta0, gp_theta1, theta0, theta1, LML


def get_pddf_acf_data(folder):
    data = np.loadtxt(f"{folder}/data_pddf_acf.txt", skiprows=1, delimiter=",", unpack=True)
    z, pddf, all_acf = data[0,:], data[1,:], data[2:,:]
    return z, pddf, all_acf


def plot_LML_contour(tex_lw=240.71031, ppi=72):
    #folder = "../data/20241125_rand"
    #folder = "../data/20241202_andes_rand"
    folder = "../data/20241206_andes_rand"
    # all_feature_name = ["Rmu", "n", "sigma", "sqrtD", "gxy"]
    # all_feature_tex = [r"$R_\mu/R_0$", r"$n$", r"$R_\sigma$", r"$\sqrt{D}/R_0$", r"$\gamma/R_0$"]
    feature_to_tex = {"sigma": r"$R_s$", "sqrtD": r"$D_2$", "gxy": r"$\gamma L$"}
    feature_to_ACF_tex = {"sigma": r"$C_{R_s}(z)$", "sqrtD": r"$C_{D_2}(z)$", "gxy": r"$C_{\gamma L}(z)$"}

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax0 = fig.add_subplot(221)

    ax1 = fig.add_subplot(222)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)

    # plot pddf and acf
    z, pddf, all_acf = get_pddf_acf_data(folder)
    all_feature_name = ["Rmu", "n", "sigma", "sqrtD", "gxy"]

    for feature_name, feature_tex in feature_to_tex.items():
        feature_index = all_feature_name.index(feature_name)
        ax0.plot(z, all_acf[feature_index], "-", lw=1, label= feature_tex)

    ax0.plot(z, pddf, "-", lw=1, label=r"$p(z)$")
    ax0.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax0.xaxis.set_minor_locator(plt.MultipleLocator(0.5))

    ax0.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax0.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax0.set_ylim(-1, 1.6)

    ax0.set_xlabel(r"$z$", fontsize=9, labelpad=0)
    ax0.set_ylabel(r"$ACF$", fontsize=9, labelpad=0)
    ax0.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax0.legend(loc="upper right", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.5, frameon=False, fontsize=9)

    axs = [ax1, ax2, ax3]

    grid_size = 2

    ticks = [
        (np.logspace(-1, 0, grid_size), np.logspace(-3, -1, grid_size)),
        (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
        (np.logspace(-1, 1, grid_size), np.logspace(-4, -2, grid_size)),
    ]

    ticklabels = [
        ((r"$10^{-1}$", r"$10^{0}$"), (r"$10^{-3}$", r"$10^{-1}$")),
        ((r"$10^{-1}$", r"$10^{0}$"), (r"$10^{-3}$", r"$10^{-2}$")),
        ((r"$10^{-1}$", r"$10$"), (r"$10^{-4}$", r"$10^{-2}$")),
    ]

    i = -1
    for feature_name, feature_tex in feature_to_tex.items():
        i += 1
        gp_theta0, gp_theta1, theta0, theta1, LML = get_LML_date(folder, feature_name)
        Theta0, Theta1 = np.meshgrid(theta0, theta1)
        axs[i].contour(Theta0, Theta1, LML, levels=50, linewidths=1)  # , cmap="summer")
        axs[i].plot([gp_theta0[0]], [gp_theta1[0]], "x", color="black", markersize=5, markeredgewidth=1)  # , label=r"l=%.2e, $\sigma$=%.2e" % (gp_theta0[0], gp_theta1[0]))
        axs[i].set_xscale("log")
        axs[i].set_yscale("log")
        axs[i].tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
        # for tick in axs[i].get_yticklabels():
        #    tick.set_rotation(90)
        axs[i].set_xlabel(r"$l$", labelpad=0, fontsize=9)
        axs[i].set_ylabel(r"$\sigma$", labelpad=0, fontsize=9)
        #axs[i].set_xticks(ticks[i][0])
        #axs[i].set_yticks(ticks[i][1])
        #axs[i].set_xticklabels(ticklabels[i][0])
        #axs[i].set_yticklabels(ticklabels[i][1])
        axs[i].minorticks_off()

        """
        yticks = axs[i].yaxis.get_minor_ticks()
        for tick in yticks:
            if np.log10(tick.get_loc())%1 != 0:
                tick.label1.set_visible(False)
        xticks = axs[i].xaxis.get_minor_ticks()
        for tick in xticks:
            if np.log10(tick.get_loc())%1 != 0:
                tick.label1.set_visible(False)
        """
        #axs[i].legend(title=feature_tex, loc="upper left", fontsize=9, frameon=False)
        #axs[i].set_title(feature_tex, fontsize=9)
        axs[i].legend(title=feature_tex, fontsize=9, frameon=False)

    x0, y0 = 0.3, 0.9
    x1, y1 = 0.2, 0.5

    # Plot arrow on ax2 using converted scale
    ax_add = fig.add_subplot(224)
    ax_add.set_xlim(0,1)
    ax_add.set_ylim(0,1)
    ax_add.arrow(x0, y0, x1 - x0, y1 - y0, head_width=0.05, head_length=0.05, fc="k", ec="k", lw=1)
    ax_add.text(0.35, 0.7, r"LML increase", fontsize=9, ha="center", va="center", rotation=80)
    ax_add.set_axis_off()

    # cbar = fig.colorbar(axs[i].collections[0], ax=axs[i], fraction=0.046, pad=0.04)
    # cbar.ax.tick_params(labelsize=7)

    # axall = fig.add_subplot(111, frameon=False)
    # axall.tick_params(labelcolor="none", which="both", top=False, bottom=False, left=False, right=False)
    # axall.set_xlabel(r"$l$", fontsize=9)
    # axall.set_ylabel(r"$\sigma$", fontsize=9)
    # ax13.set_ylabel("ML Inversion", fontsize=9, labelpad=0)

    ax0.text(0.1, 0.1, r"$(a)$", fontsize=9, transform=ax0.transAxes)
    ax1.text(0.1, 0.1, r"$(b)$", fontsize=9, transform=ax1.transAxes)
    ax2.text(0.1, 0.1, r"$(c)$", fontsize=9, transform=ax2.transAxes)
    ax3.text(0.1, 0.1, r"$(d)$", fontsize=9, transform=ax3.transAxes)

    plt.tight_layout(pad=0.3)
    # plt.subplots_adjust(left=0.12, bottom=0.1)

    plt.savefig("figures/LML_contour.pdf", format="pdf", dpi=300)
    plt.savefig("figures/LML_contour.png", dpi=300)
    plt.show()
    plt.close()


def get_GRP_data(folder, feature):
    data = np.loadtxt(f"{folder}/data_{feature}_prediction.txt", skiprows=1, delimiter=",", unpack=True)
    return data


def plot_GPR_prediction(tex_lw=240.71031, ppi=72):
    #folder = "../data/20241125_rand"
    folder = "../data/20241206_andes_rand"

    all_feature, sqv0, sqv1, sqv2 = read_SVD_data(folder)
    feature_to_tex = {"sigma": r"$R_s$", "sqrtD": r"$D_2$", "gxy": r"$\gamma L$"}

    fig = plt.figure(figsize=(tex_lw / ppi * 2, tex_lw / ppi * 0.7))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    axs = [ax1, ax2, ax3]

    major_locator = [0.1, 1, 10]
    minor_locator = [0.05, 0.5, 5]
    i = -1
    for feature_name, feature_tex in feature_to_tex.items():
        i += 1
        mu, mu_pred, mu_err = get_GRP_data(folder, feature_name)

        # axs[i].scatter(mu, mu_pred, s=0.5, marker=".", c=relative_err, cmap="rainbow", norm=norm)
        axs[i].scatter(mu, mu_pred, s=0.5, marker=".", c="royalblue")
        err = np.average(np.abs(mu - mu_pred) / np.maximum(mu, mu_pred))
        axs[i].text(0.6, 0.3, f"Err={err:.2f}", transform=axs[i].transAxes, fontsize=9)
        axs[i].plot(mu, mu, color="gray", linestyle="--", lw=0.25, alpha=0.5)
        axs[i].xaxis.set_major_locator(plt.MultipleLocator(major_locator[i]))
        axs[i].xaxis.set_minor_locator(plt.MultipleLocator(minor_locator[i]))
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(major_locator[i]))
        axs[i].yaxis.set_minor_locator(plt.MultipleLocator(minor_locator[i]))
        axs[i].grid(True, which="major", linestyle="--", linewidth=0.5)
        axs[i].tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
        # axs[i].legend(title = feature_tex, fontsize=9, loc="upper left")
        axs[i].text(0.2, 0.6, feature_tex, transform=axs[i].transAxes, fontsize=9)

    axall = fig.add_subplot(111, frameon=False)
    axall.tick_params(labelcolor="none", which="both", top=False, bottom=False, left=False, right=False)
    axall.set_xlabel("MC References", fontsize=9, labelpad=-3)
    axall.set_ylabel("ML Inversion", fontsize=9, labelpad=-3)
    # ax13.set_ylabel("ML Inversion", fontsize=9, labelpad=0)

    ax1.text(0.8, 0.1, r"$(a)$", fontsize=9, transform=ax1.transAxes)
    ax2.text(0.8, 0.1, r"$(b)$", fontsize=9, transform=ax2.transAxes)
    ax3.text(0.8, 0.1, r"$(c)$", fontsize=9, transform=ax3.transAxes)

    plt.tight_layout(pad=0.5)
    # plt.subplots_adjust(left=0.08, bottom=0.2)
    plt.savefig("figures/GPR_prediction.pdf", format="pdf", dpi=300)
    plt.savefig("figures/GPR_prediction.png", dpi=300)
    plt.show()
    plt.close()
