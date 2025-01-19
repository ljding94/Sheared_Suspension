import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
import os


def get_feature_Iq_gq_data(folder, finfos, save_IqIq_af=False):
    all_Rmu = []  # system size related
    all_n = []  # system size related
    all_sigma = []
    all_sqrtD = []
    all_gxy = []

    all_Iq2D = []
    all_Iq2D_af = []
    all_IqIq_af = []
    all_gq = []
    all_gq_phi = []
    all_gq_r = []

    qr = []
    qphi = []
    for i in range(len(finfos)):
        filename = f"{folder}/obs_{finfos[i]}.csv"

        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        bnum_phi = len(data[0]) - 7
        print(f"bnum_phi: {bnum_phi}")
        bnum_r = int((len(data) - 2) / 2)
        if save_IqIq_af:
            bnum_r = int((len(data) - 2) / 3)
        print(f"bnum_r: {bnum_r}")

        Rmu, n, sigma, sqrtD, gxy = data[0, 1], data[0, 2], data[0, 3], data[0, 4], data[0, 5]

        qphi = data[1, 7:]
        R0 = 0.00125
        qr = R0 * data[2 : 2 + bnum_r, 6]
        Iq2D = data[2 : 2 + bnum_r, 7:]
        IqIq_af = data[2 + bnum_r : 2 + 2 * bnum_r, 7:]

        gq = IqIq_af / (Iq2D * Iq2D)
        QPHI, QR = np.meshgrid(qphi, qr)
        # gq_phi = np.mean(gq*QR, axis=0)/np.mean(QR, axis=0)
        gq_phi = np.mean(gq, axis=0)
        gq_r = np.mean(gq, axis=1)

        # q phi is from [0,pi], fill [0,2pi], and
        # qphi = np.concatenate((qphi, qphi[1:] + np.pi))
        # also fill Iq2D using symmetry
        Iq2D = np.concatenate((Iq2D, Iq2D[:, 1:]), axis=1)
        IqIq_af = np.concatenate((IqIq_af, IqIq_af[:, 1:]), axis=1)

        if save_IqIq_af:
            Iq2D_af = data[2 + 2 * bnum_r : 2 + 3 * bnum_r, 7:]
            Iq2D_af = np.concatenate((Iq2D_af, Iq2D_af[:, 1:]), axis=1)
            all_Iq2D_af.append(Iq2D_af)

        gq = np.concatenate((gq, gq[:, 1:]), axis=1)

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
    all_feature_tex = [r"$R_\mu/R_0$", r"$nL^2$", r"$R_s$", r"$D_2$", r"$\gamma L$"]

    qphi = np.array(qphi)
    qr = np.array(qr)
    all_Iq2D = np.array(all_Iq2D)
    all_Iq2D_af = np.array(all_Iq2D_af)
    all_IqIq_af = np.array(all_IqIq_af)
    all_gq = np.array(all_gq)
    all_gq_phi = np.array(all_gq_phi)
    all_gq_r = np.array(all_gq_r)

    if save_IqIq_af:
        return all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_Iq2D_af, all_IqIq_af, all_gq, all_gq_r, all_gq_phi, qr, qphi
    else:
        return all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi


def plot_config_Iq(tex_lw=240.71031, ppi=72):
    folder = "../data/data_local/20241202"
    finfo = "n151_Rmu1.0_sigma0.0_sqrtD0.5_gxy30.0"

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1.25))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    # ax0 = plt.subplot(141)
    # ax1 = plt.subplot(142, projection="polar")
    # ax2 = plt.subplot(143, projection="polar")
    # ax3 = plt.subplot(144, projection="polar")

    ax0 = plt.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
    ax1 = plt.subplot2grid((3, 3), (2, 0), projection="polar")
    ax2 = plt.subplot2grid((3, 3), (2, 1), projection="polar")
    ax3 = plt.subplot2grid((3, 3), (2, 2), projection="polar")

    # plot config

    filename = f"{folder}/config_{finfo}.csv"
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)
    R0, R, x, y, x_af, y_af = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]
    # Plot the data
    L = 1.0 / R0[0]
    ax0.plot(np.array([-0.5, 0.5, 0.5, -0.5, -0.5]) * L, np.array([-0.5, -0.5, 0.5, 0.5, -0.5]) * L, color="black", lw=2)  # frame

    s = 200 * (1 * ax0.get_window_extent().width / (1.1 * L + 1.0) * 72.0 / fig.dpi) ** 2
    ax0.scatter(x / R0, y / R0, s=s * R / R0, facecolors="royalblue", edgecolors="none", label=r"$\mathcal{S}$")
    ax0.scatter(x_af / R0, y_af / R0, s=s * R / R0, facecolors="tomato", edgecolors="none", label=r"$\Gamma(\mathcal{S})$")
    ax0.set_xlim(-0.6 * L, 0.6 * L)
    ax0.set_ylim(-0.6 * L, 0.6 * L)

    ax0.plot([0, 0], [0, 0.45 * L], "-", color="gray", lw=1, zorder=-1)
    ax0.plot([0, 0.225 * L], [0, 0.45 * L], "-", color="gray", lw=1, zorder=-1)
    ax0.arrow(0, 0.2 * L, 0.1 * L - 0.03 * L, 0, head_width=0.02 * L, head_length=0.02 * L, fc="black", ec="black", lw=1)
    ax0.arrow(0, 0.3 * L, 0.15 * L - 0.03 * L, 0, head_width=0.02 * L, head_length=0.02 * L, fc="black", ec="black", lw=1)
    ax0.arrow(0, 0.4 * L, 0.2 * L - 0.03 * L, 0, head_width=0.02 * L, head_length=0.02 * L, fc="black", ec="black", lw=1)

    ax0.set_aspect("equal")
    ax0.legend(loc="center right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, fontsize=9)
    ax0.set_axis_off()

    # Iq
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_Iq2D_af, all_IqIq_af, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, [finfo], True)
    qphi = np.concatenate((qphi, qphi[1:] + np.pi))
    QPHI, QR = np.meshgrid(qphi, qr)

    p1 = ax1.pcolormesh(QPHI, QR, np.log10(all_Iq2D[0]), cmap="rainbow", shading="gouraud", rasterized=True)
    cbar1 = plt.colorbar(p1, ax=ax1, orientation="vertical", shrink=0.5, pad=0.02)
    cbar1.ax.xaxis.set_ticks_position("top")
    cbar1.ax.xaxis.set_label_position("top")
    cbar1.ax.set_title(r"$log_{10}I(\vb{q,\mathcal{S}})$", fontsize=9, x=-7)
    cbar1.ax.tick_params(direction="in", labelsize=7)
    cbar1.ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax1.set_axis_off()

    p2 = ax2.pcolormesh(QPHI, QR, np.log10(all_Iq2D_af[0]), cmap="rainbow", shading="gouraud", rasterized=True)
    cbar2 = plt.colorbar(p2, ax=ax2, orientation="vertical", shrink=0.5, pad=0.02)
    cbar2.ax.xaxis.set_ticks_position("top")
    cbar2.ax.xaxis.set_label_position("top")
    cbar2.ax.set_title(r"$log_{10}I(\vb{q,\Gamma(\mathcal{S})})$", fontsize=9, x=-7)
    cbar2.ax.tick_params(direction="in", labelsize=7, pad=1.2)
    cbar2.ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax2.set_axis_off()

    #p3 = ax3.pcolormesh(QPHI, QR, 4*all_IqIq_af[0]/(all_Iq2D_af[0]+all_Iq2D[0])**2, cmap="rainbow", shading="gouraud")
    #p3 = ax3.pcolormesh(QPHI, QR, all_IqIq_af[0]/all_Iq2D[0]**2, cmap="rainbow", shading="gouraud")
    p3 = ax3.pcolormesh(QPHI, QR, all_IqIq_af[0], cmap="rainbow", shading="gouraud", vmax=5, rasterized=True)
    # p3 = ax3.pcolormesh(QPHI, QR, all_IqIq_af[0]/all_Iq2D[0], cmap="rainbow", shading="gouraud", vmax=5, rasterized=True)
    cbar3 = plt.colorbar(p3, ax=ax3, orientation="vertical", shrink=0.5, pad=0.02)
    cbar3.ax.xaxis.set_ticks_position("top")
    cbar3.ax.xaxis.set_label_position("top")
    #cbar3.ax.set_title(r"$I(\vb{q,\Gamma(\mathcal{S})})I(\vb{q,\mathcal{S}})$", fontsize=9, x=-7, y=1.25)
    cbar3.ax.set_title(r"$I(\vb{q,\Gamma(\mathcal{S})})I(\vb{q,\mathcal{S}})$", fontsize=9, x=-7)
    cbar3.ax.tick_params(direction="in", labelsize=7, pad=1.2)
    # cbar3.ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax3.set_axis_off()
    # ax3.set_title(r"$I(\vb{q,\Gamma(\mathcal{S})})I(\vb{q,\mathcal{S}})$", fontsize=9, y=0.95)

    ax0.text(0.95, 0.2, r"$(a)$", fontsize=9, transform=ax0.transAxes)
    ax1.text(0.8, 0.0, r"$(b)$", fontsize=9, transform=ax1.transAxes)
    ax2.text(0.8, 0.0, r"$(c)$", fontsize=9, transform=ax2.transAxes)
    ax3.text(0.8, 0.0, r"$(d)$", fontsize=9, transform=ax3.transAxes)

    plt.tight_layout(pad=0.05)

    plt.savefig("./figures/config_Iq.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/config_Iq.png", format="png", dpi=300)
    plt.show()
    plt.close()


def plot_Iq_gq_single(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.35))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((1, 3), (0, 0), projection="polar")
    ax01 = plt.subplot2grid((1, 3), (0, 1), sharex=ax00, projection="polar")
    ax02 = plt.subplot2grid((1, 3), (0, 2), sharex=ax00, projection="polar")

    # folder = "../data/data_local/20241125"
    # folder = "../data/20241129_prec"
    # finfo = "n100_Rmu1.0_sigma0.0_sqrtD1.0_gxy10.0"

    folder = "../data/20241202_andes_prec"
    finfo = "n150_Rmu1.0_sigma0.0_sqrtD0.5_gxy10.0"


    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, [finfo])
    # print(qphi)
    qphi = np.concatenate((qphi, qphi[1:] + np.pi))
    # print(qphi)

    Iq2D = all_Iq2D[0]
    gq = all_gq[0]
    QPHI, QR = np.meshgrid(qphi, qr)

    gq_theory = calc_gq(QR, QPHI, 0.5, 10.0)

    p00 = ax00.pcolormesh(QPHI, QR, np.log10(Iq2D), cmap="rainbow", shading="gouraud", rasterized=True)
    # ax00.set_yscale("log")
    # Cs = ax00.contour(QPHI, QR, np.log10(Iq2D), colors="gray", linewidths=0.5, linestyle=":")
    # ax00.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.1f", colors="black")
    # cbar1 = plt.colorbar(p1, ax=ax1, orientation="vertical", shrink=0.5)
    #cbar00 = plt.colorbar(p00, ax=ax00, orientation="vertical", shrink=0.5, pad=0.02)
    cbar00 = plt.colorbar(p00, ax=ax00, orientation="vertical", shrink=0.5, pad=0.02)
    cbar00.ax.xaxis.set_ticks_position("top")
    cbar00.ax.xaxis.set_label_position("top")
    cbar00.ax.set_title(r"$log_{10}I(\vb{q})$", fontsize=9)
    cbar00.ax.tick_params(direction="in", labelsize=7, pad=1.2)
    cbar00.ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    # cbar1.set_label(r"$log_{10}(\left<I(q)\right>)$")

    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    # ax00.set_yticks([1000, 3000])
    # ax00.set_yticklabels([ r"$1e3$", r"$3e3$"], fontsize=7)
    ax00.grid(axis="y", linewidth=0.5, linestyle="--")
    ax00.grid(axis="x", linestyle="None")
    # ax00.set_rlabel_position(160)
    # ax00.set_aspect("equal")
    ax00.set_axis_off()
    # ax00.set_ylabel(r"$\theta$", fontsize=9)
    # ax00.set_xlabel(r"$q$", fontsize=9)
    # ax00.set_title(r"$log_{10}(\left<I(\vb{q})\right>)$", fontsize=9)

    # p3 = ax3.pcolormesh(QPHI, QR, np.log10(IqIq_af), cmap="rainbow", shading="gouraud")
    p01 = ax01.pcolormesh(QPHI, QR, gq, cmap="rainbow", shading="gouraud", rasterized=True, vmin=1, vmax=1.9)
    # ax01.set_yscale("log")
    # Cs = ax01.contour(QPHI, QR, IqIq_af, colors="gray", linewidths=0.5, linestyle=":")
    # ax01.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.3f", colors="black")
    # cbar01 = plt.colorbar(p01, ax=ax01, orientation="horizontal", location='top', shrink=0.5, pad=0.0)
    cbar01 = plt.colorbar(p01, ax=ax01, orientation="vertical", shrink=0.5, pad=0.02)
    cbar01.ax.xaxis.set_ticks_position("top")
    cbar01.ax.xaxis.set_label_position("top")
    cbar01.ax.set_title(r"$g(\vb{q})$", fontsize=9)
    cbar01.ax.tick_params(direction="in", labelsize=7, pad=1.2)
    cbar01.ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    cbar01.ax.set_ylim(1,1.9)
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)

    # ax01.set_yticks([1000, 3000])
    # ax01.set_yticklabels([r"$1e3$", r"$3e3$"], fontsize=7)
    # ax01.set_rlabel_position(160)
    ax01.grid(axis="y", linewidth=0.5, linestyle="--")
    ax01.grid(axis="x", linestyle="None")

    # ax01.set_xlabel(r"$\theta$", fontsize=9)
    # ax01.set_ylabel(r"$q$", fontsize=9)
    ax01.set_axis_off()


    p02 = ax02.pcolormesh(QPHI, QR, gq_theory, cmap="rainbow", shading="gouraud", rasterized=True, vmin=1, vmax=1.9)
    # ax02.set_yscale("log")
    # Cs = ax02.contour(QPHI, QR, IqIq_af, colors="gray", linewidths=0.5, linestyle=":")
    # ax02.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.3f", colors="black")
    # cbar02 = plt.colorbar(p02, ax=ax02, orientation="horizontal", location='top', shrink=0.5, pad=0.0)
    cbar02 = plt.colorbar(p02, ax=ax02, orientation="vertical", shrink=0.5, pad=0.02)
    cbar02.ax.xaxis.set_ticks_position("top")
    cbar02.ax.xaxis.set_label_position("top")
    cbar02.ax.set_title(r"$g(\vb{q})$", fontsize=9)
    cbar02.ax.tick_params(direction="in", labelsize=7, pad=1.2)
    cbar02.ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    cbar02.ax.set_ylim(1, 1.9)
    ax02.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)

    # ax02.set_yticks([1000, 3000])
    # ax02.set_yticklabels([r"$1e3$", r"$3e3$"], fontsize=7)
    # ax02.set_rlabel_position(160)
    ax02.grid(axis="y", linewidth=0.5, linestyle="--")
    ax02.grid(axis="x", linestyle="None")

    # ax02.set_xlabel(r"$\theta$", fontsize=9)
    # ax02.set_ylabel(r"$q$", fontsize=9)
    ax02.set_axis_off()
    ax02.set_title("theory", fontsize=9, y=0.92)


    ax00.text(0.85, 0.05, r"$(a)$", fontsize=9, transform=ax00.transAxes)
    ax01.text(0.85, 0.05, r"$(b)$", fontsize=9, transform=ax01.transAxes)
    ax02.text(0.85, 0.05, r"$(c)$", fontsize=9, transform=ax02.transAxes)

    for ax in [ax00, ax01]:
        pass
        # ax.set_thetamin(0)
        # ax.set_thetamax(180)
        # ax.grid(False)

        # ax.set_axis_off()
    plt.tight_layout(pad=0.15)
    plt.savefig("./figures/Iq_gq_single.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/Iq_gq_single.png", format="png", dpi=300)
    plt.show()




def plot_Iq_gq_avg(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1.2))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((3, 2), (0, 0))
    ax01 = plt.subplot2grid((3, 2), (0, 1))

    ax10 = plt.subplot2grid((3, 2), (1, 0))
    ax11 = plt.subplot2grid((3, 2), (1, 1))

    ax20 = plt.subplot2grid((3, 2), (2, 0))
    ax21 = plt.subplot2grid((3, 2), (2, 1))

    # folder = "../data/data_local/20241125"
    # folder = "../data/20241130_prec"

    folder = "../data/20241202_andes_prec"
    finfo = "n150_Rmu1.0_sigma0.0_sqrtD1.0_gxy10.0"

    gq_r_tex = r"$g(q)$"
    gq_phi_tex = r"$g(\theta)$"
    # q_ticks = [500, 1000, 2000, 3000]
    q_ticks = [1, 2, 3, 5]
    q_ticks_label = [r"$1$", r"$2$", r"$3$", r"$5$"]
    # q_ticks_label = [r"$5\times 10^2$", r"$1\times 10^3$", r"$2\times 10^3$", r"$3\times 10^3$"]
    # q_ticks_label = [r"$5\times 10^2$", r"$1\times 10^3$", r"$2\times 10^3$", r"$3\times 10^3$"]
    # q_ticks_label = [r"$5e2$", r"$1e3$", r"$2e3$", r"$3e3$"]
    phi_ticks = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
    phi_ticks_label = [r"$0$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"]

    colors = ["royalblue", "tomato", "forestgreen", "darkorange"]

    # versus sigma
    sigmas = [0.0, 0.2, 0.3]
    for i in range(len(sigmas)):
        finfo = f"n150_Rmu1.0_sigma{sigmas[i]:.1f}_sqrtD1.0_gxy10.0"

        all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, [finfo])
        gq_r = all_gq_r[0]
        gq_phi = all_gq_phi[0]

        color = colors[i]
        if(sigmas[i] == 0.0):
            color = "black"

        ax00.semilogx(qr, gq_r, lw=1, label=f"{sigmas[i]:.1f}", color=color)
        # ax00.set_xlabel(r"$q$", fontsize=9)
        ax00.set_ylabel(gq_r_tex, fontsize=9)
        ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
        ax00.set_xticks(q_ticks)
        # ax00.set_xticklabels(q_ticks_label)
        ax00.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax00.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

        ax01.plot(qphi, gq_phi, lw=1, label=f"{sigmas[i]:.1f}", color=color)
        # ax01.set_xlabel(r"$\theta$", fontsize=9)
        ax01.set_ylabel(gq_phi_tex, fontsize=9)
        ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
        ax01.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax01.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax01.set_xticks(phi_ticks)
        ax01.set_xticklabels(phi_ticks_label)

    ax00.legend(title=r"$R_s$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax01.legend(title=r"$R_s$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    # versus sqrtD
    sqrtDs = [0.5, 1.0, 2.0, 3.0]
    for i in range(len(sqrtDs)):
        finfo = f"n150_Rmu1.0_sigma0.0_sqrtD{sqrtDs[i]:.1f}_gxy10.0"
        all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, [finfo])
        gq_r = all_gq_r[0]
        gq_phi = all_gq_phi[0]

        color = colors[i]
        if(sqrtDs[i] == 1.0):
            color = "black"
        ax10.semilogx(qr, gq_r, lw=1, label=f"{sqrtDs[i]:.1f}", color=color)
        # ax10.set_xlabel(r"$q$", fontsize=9)
        ax10.set_ylabel(gq_r_tex, fontsize=9)
        ax10.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
        ax10.set_xticks(q_ticks)
        ax10.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax10.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

        ax11.plot(qphi, gq_phi, lw=1, label=f"{sqrtDs[i]:.1f}", color=color)
        # ax11.set_xlabel(r"$\theta$", fontsize=9)
        ax11.set_ylabel(gq_phi_tex, fontsize=9)
        ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
        ax11.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax11.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax11.set_xticks(phi_ticks)
        ax11.set_xticklabels(phi_ticks_label)

    ax10.legend(title=r"$D_2$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax11.legend(title=r"$D_2$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    # versus gxy
    gxys = [10.0, 20.0, 30.0]
    for i in range(len(gxys)):
        finfo = f"n150_Rmu1.0_sigma0.0_sqrtD1.0_gxy{gxys[i]:.1f}"
        all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, [finfo])
        gq_r = all_gq_r[0]
        gq_phi = all_gq_phi[0]

        color = colors[i]
        if(gxys[i] == 10.0):
            color = "black"

        ax20.semilogx(qr, gq_r, lw=1, label=f"{gxys[i]:.0f}", color=color)
        ax20.set_xlabel(r"$q$", fontsize=9)
        ax20.set_ylabel(gq_r_tex, fontsize=9)
        ax20.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
        ax20.set_xticks(q_ticks)
        ax20.set_xticklabels(q_ticks_label)
        ax20.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax20.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

        ax21.plot(qphi, gq_phi, lw=1, label=f"{gxys[i]:.0f}", color=color)
        ax21.set_xlabel(r"$\theta$", fontsize=9)
        ax21.set_ylabel(gq_phi_tex, fontsize=9)
        ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
        ax21.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax21.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
        ax21.set_xticks(phi_ticks)
        ax21.set_xticklabels(phi_ticks_label)

    ax20.legend(title=r"$\gamma L$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.legend(title=r"$\gamma L$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

    ax00.text(0.2, 0.8, r"$(a)$", fontsize=9, transform=ax00.transAxes)
    ax01.text(0.2, 0.8, r"$(b)$", fontsize=9, transform=ax01.transAxes)
    ax10.text(0.2, 0.8, r"$(c)$", fontsize=9, transform=ax10.transAxes)
    ax11.text(0.2, 0.8, r"$(d)$", fontsize=9, transform=ax11.transAxes)
    ax20.text(0.2, 0.8, r"$(e)$", fontsize=9, transform=ax20.transAxes)
    ax21.text(0.2, 0.8, r"$(f)$", fontsize=9, transform=ax21.transAxes)

    plt.tight_layout(pad=0.2)
    plt.savefig("./figures/gq_avg.pdf", format="pdf")
    plt.savefig("./figures/gq_avg.png", format="png")
    plt.show()


def calc_gq(qr, qphi, sqrtD, gamma):
    return 1 + np.exp(-qr * qr * sqrtD * sqrtD) * np.sinc(qr * gamma * np.cos(qphi) * 0.5 / np.pi) ** 2

    # note, for numpy sinc(x) = sin(pi x)/(pi x), but for physics sinc(x) = sin(x)/x

def plot_theory_gq(tex_lw=240.71031, ppi=72):

    fig = plt.figure(figsize=(tex_lw / ppi * 2, tex_lw / ppi * 1))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((2, 3), (0, 0), projection="polar")
    ax01 = plt.subplot2grid((2, 3), (0, 1), sharex=ax00, projection="polar")
    ax02 = plt.subplot2grid((2, 3), (0, 2))
    ax10 = plt.subplot2grid((2, 3), (1, 0), projection="polar")
    ax11 = plt.subplot2grid((2, 3), (1, 1), sharex=ax10, projection="polar")
    ax12 = plt.subplot2grid((2, 3), (1, 2))

    #folder = "../data/20241202_andes_prec"
    # folder = "../data/20241129_prec"
    folder = "../data/data_local/20241202"


    # only diffusion
    finfo = "n150_Rmu1.0_sigma0.0_sqrtD0.5_gxy0.0"
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, [finfo])
    # print(qphi)
    qphi_c = np.concatenate((qphi, qphi[1:] + np.pi))
    # print(qphi)

    Iq2D = all_Iq2D[0]
    gq = all_gq[0]
    QPHI, QR = np.meshgrid(qphi_c, qr)

    gq_theory = calc_gq(QR, QPHI, 0.5, 0.0)
    p00 = ax00.pcolormesh(QPHI, QR, gq_theory, cmap="rainbow", shading="gouraud", rasterized=True)
    cbar00 = plt.colorbar(p00, ax=ax00, orientation="vertical", shrink=0.5, pad=0.1)
    cbar00.ax.xaxis.set_ticks_position("top")
    cbar00.ax.xaxis.set_label_position("top")
    cbar00.ax.set_title(r"$g(\vb{q})$", fontsize=9)
    cbar00.ax.tick_params(direction="in", labelsize=7)
    cbar00.ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    ax00.set_title("theory")
    ax00.set_axis_off()

    # p3 = ax3.pcolormesh(QPHI, QR, np.log10(IqIq_af), cmap="rainbow", shading="gouraud")
    p01 = ax01.pcolormesh(QPHI, QR, gq, cmap="rainbow", shading="gouraud", rasterized=True)
    # ax01.set_yscale("log")
    # Cs = ax01.contour(QPHI, QR, IqIq_af, colors="gray", linewidths=0.5, linestyle=":")
    # ax01.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.3f", colors="black")
    # cbar01 = plt.colorbar(p01, ax=ax01, orientation="horizontal", location='top', shrink=0.5, pad=0.0)
    cbar01 = plt.colorbar(p01, ax=ax01, orientation="vertical", shrink=0.5, pad=0.1)
    cbar01.ax.xaxis.set_ticks_position("top")
    cbar01.ax.xaxis.set_label_position("top")
    cbar01.ax.set_title(r"$g(\vb{q})$", fontsize=9)
    cbar01.ax.tick_params(direction="in", labelsize=7)
    cbar01.ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    ax01.set_title("simulation")
    ax01.set_axis_off()

    gq_theory_r = np.mean(gq_theory, axis=1)
    print("gq_theory_r.shape", gq_theory_r.shape)
    ax02.semilogx(qr, gq_theory_r, lw=1, label=r"$1+\text{sinc}(\frac{q\cos(\theta) \gamma L}{6} )e^{-q^2D^2_2}$")
    ax02.semilogx(qr[::5], all_gq_r[0][::5], "o", mfc="none", lw=1, label="simulation")

    ax02.set_xlabel(r"$q$", fontsize=9)
    ax02.set_ylabel(r"$g(q)$", fontsize=9)
    ax02.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax02.legend(loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.5, frameon=False, fontsize=9)

    # only shear
    finfo = "n150_Rmu1.0_sigma0.0_sqrtD0.0_gxy10.0"
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_gq, all_gq_r, all_gq_phi, qr, qphi = get_feature_Iq_gq_data(folder, [finfo])
    # print(qphi)
    qphi_c = np.concatenate((qphi, qphi[1:] + np.pi))
    # print(qphi)

    Iq2D = all_Iq2D[0]
    gq = all_gq[0]
    QPHI, QR = np.meshgrid(qphi_c, qr)

    gq_theory = calc_gq(QR, QPHI, 0.0, 10.0/3)
    p10 = ax10.pcolormesh(QPHI, QR, gq_theory, cmap="rainbow", shading="gouraud", rasterized=True)
    cbar10 = plt.colorbar(p10, ax=ax10, orientation="vertical", shrink=0.5, pad=0.1)
    cbar10.ax.xaxis.set_ticks_position("top")
    cbar10.ax.xaxis.set_label_position("top")
    cbar10.ax.set_title(r"$g(\vb{q})$", fontsize=9)
    cbar10.ax.tick_params(direction="in", labelsize=7)
    cbar10.ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax10.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    ax10.set_title("theory")
    ax10.set_axis_off()

    # p3 = ax3.pcolormesh(QPHI, QR, np.log10(IqIq_af), cmap="rainbow", shading="gouraud")
    p11 = ax11.pcolormesh(QPHI, QR, gq, cmap="rainbow", shading="gouraud", rasterized=True)
    cbar11 = plt.colorbar(p11, ax=ax11, orientation="vertical", shrink=0.5, pad=0.1)
    cbar11.ax.xaxis.set_ticks_position("top")
    cbar11.ax.xaxis.set_label_position("top")
    cbar11.ax.set_title(r"$g(\vb{q})$", fontsize=9)
    cbar11.ax.tick_params(direction="in", labelsize=7)
    cbar11.ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    ax11.set_title("simulation")
    ax11.set_axis_off()

    gq_theory_phi = np.mean(gq_theory, axis=0)
    print("gq_theory_phi.shape", gq_theory_phi.shape)
    ax12.plot(qphi, gq_theory_phi[:len(qphi)], lw=1, label=r"$1+\text{sinc}(\frac{q\cos(\theta) \gamma L}{6} )e^{-q^2D^2_2}$")
    ax12.plot(qphi[::5], all_gq_phi[0][::5], "o", mfc="none", lw=1, label="simulation")

    ax12.set_xlabel(r"$\theta$", fontsize=9)
    ax12.set_ylabel(r"$g(\theta)$", fontsize=9)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax12.legend(loc="upper right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.5, frameon=False, fontsize=9)



    plt.tight_layout(pad=0.2)

    plt.savefig("./figures/theory_gq.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/theory_gq.png", format="png", dpi=300)
    plt.show()
