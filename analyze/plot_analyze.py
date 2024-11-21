import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def plot_gas_config(filename, finfo, show=False):
    # Load the data
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)
    R0, R, x, y, x_af, y_af = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]
    # Plot the data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    L = 1.0 / R0[0]
    ax.plot(np.array([-0.5, 0.5, 0.5, -0.5, -0.5]) * L, np.array([-0.5, -0.5, 0.5, 0.5, -0.5]) * L, color="black", lw=0.1)  # frame

    s = (1 * ax.get_window_extent().width / (1.1 * L + 1.0) * 72.0 / fig.dpi) ** 2
    ax.scatter(x / R0, y / R0, s=s * R / R0, facecolors="royalblue", edgecolors="none", label="initial")
    ax.scatter(x_af / R0, y_af / R0, s=s * R / R0, facecolors="tomato", edgecolors="none", label="affine-transformed")
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.set_aspect("equal")
    # ax.legend()
    ax.set_axis_off()

    ax.set_title(finfo)
    plt.savefig(filename.replace(".csv", ".pdf"), format="pdf")
    if show:
        plt.show()
    plt.close()


def get_feature_Sq2D_SqSq2D_data(folder, parameters, random=False):

    all_N = []  # system size related
    all_theta, all_Sx, all_phi = [], [], []  # affine transformation parameters
    all_Sq2D = []
    all_Sq2D_af = []
    all_SqSq2D = []
    all_finfo = []
    qD = []
    for i in range(len(parameters)):
        if random:
            print("to be developed")
        else:
            N, sigma, theta, Sx, phi = parameters[i]
            finfo = f"N{N:.0f}_sigma{sigma:.1f}_theta{theta:.2f}_Sx{Sx:.2f}_phi{phi:.2f}"
            filename = f"{folder}/obs_{finfo}.csv"

        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        bin_num = len(data[0]) - 10
        print(f"bin_num: {bin_num}")

        qD = data[1, 10:]
        Sq2D = data[2 : bin_num + 2, 10:]
        Sq2D_af = data[bin_num + 2 : 2 * bin_num + 2, 10:]
        SqSq2D = data[2 * bin_num + 2 : 3 * bin_num + 2, 10:]

        # Apply mask to Sq2D
        center = bin_num // 2
        mask_range = 5
        Sq2D[center - mask_range : center + mask_range + 1, center - mask_range : center + mask_range + 1] = np.nan
        Sq2D_af[center - mask_range : center + mask_range + 1, center - mask_range : center + mask_range + 1] = np.nan
        SqSq2D[center - mask_range : center + mask_range + 1, center - mask_range : center + mask_range + 1] = np.nan

        all_Sq2D.append(Sq2D)
        all_Sq2D_af.append(Sq2D_af)
        all_SqSq2D.append(SqSq2D)

        all_N.append(N)
        all_theta.append(theta)
        all_Sx.append(Sx)
        all_phi.append(phi)
        all_finfo.append(finfo)

    all_feature = np.array([all_N, all_theta, all_Sx, all_phi]).T
    all_feature_name = ["N", "theta", "Sx", "phi"]
    qD = np.array(qD)
    all_Sq2D = np.array(all_Sq2D)
    all_Sq2D_af = np.array(all_Sq2D_af)
    all_SqSq2D = np.array(all_SqSq2D)

    return all_finfo, all_feature, all_feature_name, all_Sq2D, all_Sq2D_af, all_SqSq2D, qD


def plot_gas_Sq_SqSq(folder, parameter, show=False):

    all_finfo, all_feature, all_feature_name, all_Sq2D, all_Sq2D_af, all_SqSq2D, qD = get_feature_Sq2D_SqSq2D_data(folder, [parameter])

    plt.figure(figsize=(9, 3))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133)

    Sq2D = all_Sq2D[0]
    Sq2D_af = all_Sq2D_af[0]
    SqSq2D = all_SqSq2D[0]

    qDx, qDy = np.meshgrid(qD, qD)
    # vmax=0, vmin=-2.5
    # ax1.pcolormesh(qDx, qDy, np.log10(Sq2D), vmax=0, vmin=-2.5, cmap="rainbow", shading="gouraud")
    ax1.pcolormesh(qDx, qDy, np.log10(Sq2D), cmap="rainbow", shading="gouraud")
    # Cs = ax1.contour(qDx, qDy, np.log10(Sq2D), vmax=0, vmin=-2.5, levels=np.linspace(-2.5, 0, 6), colors="gray", linewidths=0.5, linestyle=":")
    Cs = ax1.contour(qDx, qDy, np.log10(Sq2D), colors="gray", linewidths=0.5, linestyle=":")
    # ax1.pcolormesh(qDx, qDy, Sq2D, cmap="rainbow", shading='gouraud')
    # Cs = ax1.contour(qDx, qDy, Sq2D, colors="gray", linewidths=0.5, linestyle=":")
    # Cs = ax1.contour(qDx, qDy, np.log10(Sq2D), vmax=0, vmin=-2.5, levels=np.linspace(-2.5, 0, 6), colors="gray", linewidths=0.5, linestyle=":")
    ax1.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.1f", colors="black")

    ax1.set_aspect("equal")
    ax1.set_xlabel(r"$q_x$")
    ax1.set_ylabel(r"$q_y$")
    ax1.set_title(r"$log_{10}(Sq2D)$")

    # ax2.pcolormesh(qDx, qDy, np.log10(Sq2D_af), vmax=0, vmin=-2.5, cmap="rainbow", shading="gouraud")
    ax2.pcolormesh(qDx, qDy, np.log10(Sq2D_af), cmap="rainbow", shading="gouraud")
    # Cs = ax2.contour(qDx, qDy, np.log10(Sq2D_af), vmax=0, vmin=-2.5, levels=np.linspace(-2.5, 0, 6), colors="gray", linewidths=0.5, linestyle=":")
    Cs = ax2.contour(qDx, qDy, np.log10(Sq2D_af), colors="gray", linewidths=0.5, linestyle=":")
    ax2.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.1f", colors="black")
    ax2.set_aspect("equal")
    ax2.set_xlabel(r"$q_x$")
    ax2.set_ylabel(r"$q_y$")
    ax2.set_title(r"$log_{10}(Sq2D_{af})$")

    # ax3.pcolormesh(qDx, qDy, np.log10(SqSq2D), vmax=0, vmin=-1, cmap="rainbow", shading="gouraud")
    # Cs = ax3.contour(qDx, qDy, np.log10(SqSq2D), vmax=0, vmin=-1, levels=np.linspace(-1, 0, 6), colors="gray", linewidths=0.5, linestyle=":")
    # ax3.pcolormesh(qDx, qDy, SqSq2D, cmap="rainbow", vmax=1, vmin=0, shading='gouraud')
    # Cs = ax3.contour(qDx, qDy, SqSq2D, colors="gray", vmax=1, vmin=0, levels=np.linspace(0, 1, 6), linewidths=0.5, linestyle=":")
    ax3.pcolormesh(qDx, qDy, SqSq2D, cmap="rainbow", shading="gouraud")
    Cs = ax3.contour(qDx, qDy, SqSq2D, colors="gray", linewidths=0.5, linestyle=":")
    ax3.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.3f", colors="black")
    ax3.set_aspect("equal")
    ax3.set_xlabel(r"$\Delta q_x$")
    ax3.set_ylabel(r"$\Delta q_y$")
    # ax3.set_title(r"$log_{10}(SqSq2D)$")
    ax3.set_title(r"$SqSq2D$")

    plt.tight_layout()
    N, sigma, theta, Sx, phi = parameter
    finfo = f"N{N:.0f}_sigma{sigma:.1f}_theta{theta:.1f}_Sx{Sx:.1f}_phi{phi:.1f}"
    plt.savefig(f"{folder}/Sq_SqSq_{finfo}.png")
    if show:
        plt.show()
    plt.close()


def get_feature_Iq2D_IqIq_af_data(folder, finfos):

    all_R0 = []  # system size related
    all_n = []  # system size related
    all_sigma = []
    all_sqrtD = []
    all_gxy = []

    all_Iq2D = []
    all_Iq2D_af = []
    all_IqIq_af = []
    all_finfo = []
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
        bnum_r = int((len(data) - 2) / 3)
        print(f"bnum_r: {bnum_r}")

        R0, n, sigma, sqrtD, gxy = data[0, 1], data[0, 2], data[0, 3], data[0, 4], data[0, 5]

        qphi = data[1, 7:]
        qr = data[2 : 2 + bnum_r, 6]
        Iq2D = data[2 : 2 + bnum_r, 7:]
        Iq2D_af = data[2 + bnum_r : 2 + 2 * bnum_r, 7:]
        IqIq_af = data[2 + 2 * bnum_r : 2 + 3 * bnum_r, 7:]

        all_Iq2D.append(Iq2D)
        all_Iq2D_af.append(Iq2D_af)
        all_IqIq_af.append(IqIq_af)

        all_R0.append(R0)
        all_n.append(n)
        all_sigma.append(sigma)
        all_sqrtD.append(sqrtD)
        all_gxy.append(gxy)

    all_feature = np.array([all_R0, all_n, all_sigma, all_sqrtD, all_gxy]).T
    all_feature_name = ["R0", "n", "sigma", "sqrtD", "gxy"]
    all_feature_tex = [r"$R_0$", r"$n$", r"$\sigma$", r"$\sqrt{D}$", r"$\gamma_{xy}$"]

    qphi = np.array(qphi)
    qr = np.array(qr)
    all_Iq2D = np.array(all_Iq2D)
    all_Iq2D_af = np.array(all_Iq2D_af)
    all_IqIq_af = np.array(all_IqIq_af)

    return all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_Iq2D_af, all_IqIq_af, qr, qphi


def plot_gas_Iq_IqIq(folder, finfo, show=True):

    # get single data
    all_feature, all_feature_name, all_feature_tex, all_Iq2D, all_Iq2D_af, all_IqIq_af, qr, qphi = get_feature_Iq2D_IqIq_af_data(folder, [finfo])

    plt.figure(figsize=(9, 3))
    ax1 = plt.subplot(231, projection="polar")
    ax2 = plt.subplot(232, projection="polar")
    ax3 = plt.subplot(233, projection="polar")

    ax4 = plt.subplot(234)
    ax5 = plt.subplot(235)

    Iq2D = all_Iq2D[0]
    Iq2D_af = all_Iq2D_af[0]
    IqIq_af = all_IqIq_af[0]

    print(f"qr.shape: {qr.shape}")
    print(f"qphi.shape: {qphi.shape}")
    print("qr", qr)
    print("qphi", qphi)

    QPHI, QR = np.meshgrid(qphi, qr)
    # QR, QPHI = np.meshgrid(qr, qphi)
    print(f"QR.shape: {QR.shape}")
    print(f"QPHI.shape: {QPHI.shape}")
    print(f"Iq2D.shape: {Iq2D.shape}")

    p1 = ax1.pcolormesh(QPHI, QR, np.log10(Iq2D), cmap="rainbow", shading="gouraud")
    # Cs = ax1.contour(QPHI, QR, np.log10(Iq2D), colors="gray", linewidths=0.5, linestyle=":")
    # ax1.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.1f", colors="black")
    cbar1 = plt.colorbar(p1, ax=ax1, orientation="vertical", shrink=0.5)
    cbar1.set_label(r"$log_{10}(\left<I(q)\right>)$")

    ax1.set_aspect("equal")
    ax1.set_xlabel(r"$q_\phi$")
    ax1.set_ylabel(r"$q_r$")
    ax1.set_title(r"$log_{10}(\left<I(q)\right>)$")

    p2 = ax2.pcolormesh(QPHI, QR, np.log10(Iq2D_af), cmap="rainbow", shading="gouraud")
    # Cs = ax2.contour(QPHI, QR, np.log10(Iq2D_af), colors="gray", linewidths=0.5, linestyle=":")
    # ax2.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.1f", colors="black")
    cbar2 = plt.colorbar(p2, ax=ax2, orientation="vertical", shrink=0.5)
    cbar2.set_label(r"$log_{10}(\left<I'(q)\right>)$")

    ax2.set_xlabel(r"$q_\phi$")
    ax2.set_ylabel(r"$q_r$")
    ax2.set_title(r"$log_{10}(\left<I'(q)\right>)$")

    gq = IqIq_af / (Iq2D * Iq2D)
    # p3 = ax3.pcolormesh(QPHI, QR, np.log10(IqIq_af), cmap="rainbow", shading="gouraud")
    p3 = ax3.pcolormesh(QPHI, QR, gq, cmap="rainbow", shading="gouraud")
    # Cs = ax3.contour(QPHI, QR, IqIq_af, colors="gray", linewidths=0.5, linestyle=":")
    # ax3.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.3f", colors="black")
    cbar3 = plt.colorbar(p3, ax=ax3, orientation="vertical", shrink=0.5)
    cbar3.set_label(r"$g(q)$")

    ax3.set_xlabel(r"$q_\phi$")
    ax3.set_ylabel(r"$q_r$")
    ax3.set_title(r"$g(q)=\left<I(q)I'(q)\right>/\left<I(q)\right>^2$")

    IqIq_af_r = np.mean(gq, axis=1)
    print(f"IqIq_af_r.shape: {IqIq_af_r.shape}")
    print(f"qr.shape: {qr.shape}")
    ax4.semilogx(qr, IqIq_af_r, color="black")
    ax4.set_xlabel(r"$q_r$")
    ax4.set_ylabel(r"$\sum_{\phi}g(q)$")

    IqIq_af_phi = np.mean(gq, axis=0)
    ax5.plot(qphi, IqIq_af_phi, color="black")
    ax5.set_xlabel(r"$q_\phi$")
    ax5.set_ylabel(r"$\sum_{r}g(q)$")

    for ax in [ax1, ax2, ax3]:
        ax.set_thetamin(0)
        ax.set_thetamax(180)
        ax.grid(False)
        ax.set_axis_off()

    ax = plt.subplot(111)
    ax.set_axis_off()
    ax.text(0.75, 0.3, f"{finfo}", fontsize=12, ha="center", va="center", transform=ax.transAxes)
    plt.tight_layout()
    plt.savefig(f"{folder}/Iq_IqIq_{finfo}.png")
    if show:
        plt.show()
    plt.close()


def calc_single_sphere_Iq(qr, R):
    return 3.0 / (qr * R) ** 3 * (np.sin(qr * R) - qr * R * np.cos(qr * R))


def test_plot():
    r = np.linspace(400, 4000, 40)
    phi = np.linspace(0, 2 * np.pi / 20 * 19, 20)
    print(phi)

    R, PHI = np.meshgrid(r, phi)

    Z = np.sin(PHI) * np.exp(-R / 1000)
    print("R.shape", R.shape)
    print("Z.shape", Z.shape)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    ax.pcolormesh(PHI, R, Z, shading="gouraud")
    ax.set_aspect("equal")
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"Polar Coordinate Plot")
    plt.show()
