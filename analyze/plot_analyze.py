import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def plot_gas_config(filename, finfo, show=False):
    # Load the data
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)
    f, x, y, x_af, y_af = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    # Plot the data
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([-0.5, 0.5, 0.5, -0.5, -0.5], [-0.5, -0.5, 0.5, 0.5, -0.5], color="black")  # frame
    ax.scatter(x, y, s=f, color="gray", label="initial")
    ax.scatter(x_af, y_af, s=f, color="tomato", label="affine-transformed")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.set_aspect("equal")
    # ax.legend()
    ax.set_axis_off()

    ax.set_title(finfo)
    plt.savefig(filename.replace(".csv", ".png"))
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
            finfo = f"N{N:.0f}_sigma{sigma:.1f}_theta{theta:.1f}_Sx{Sx:.1f}_phi{phi:.1f}"
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
    ax1.pcolormesh(qDx, qDy, np.log10(Sq2D), vmax=0, vmin=-2.5, cmap="rainbow", shading="gouraud")
    Cs = ax1.contour(qDx, qDy, np.log10(Sq2D), vmax=0, vmin=-2.5, levels=np.linspace(-2.5, 0, 6), colors="gray", linewidths=0.5, linestyle=":")
    # ax1.pcolormesh(qDx, qDy, Sq2D, cmap="rainbow", shading='gouraud')
    # Cs = ax1.contour(qDx, qDy, Sq2D, colors="gray", linewidths=0.5, linestyle=":")
    # Cs = ax1.contour(qDx, qDy, np.log10(Sq2D), vmax=0, vmin=-2.5, levels=np.linspace(-2.5, 0, 6), colors="gray", linewidths=0.5, linestyle=":")
    ax1.clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.1f", colors="black")
    ax1.set_aspect("equal")
    ax1.set_xlabel(r"$q_x$")
    ax1.set_ylabel(r"$q_y$")
    ax1.set_title(r"$log_{10}(Sq2D)$")

    ax2.pcolormesh(qDx, qDy, np.log10(Sq2D_af), vmax=0, vmin=-2.5, cmap="rainbow", shading="gouraud")
    Cs = ax2.contour(qDx, qDy, np.log10(Sq2D_af), vmax=0, vmin=-2.5, levels=np.linspace(-2.5, 0, 6), colors="gray", linewidths=0.5, linestyle=":")
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
