import json
from subprocess import PIPE, Popen

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, MultipleLocator
from mpl_toolkits.axes_grid1 import host_subplot


def show_psnr_ssim_plot(filename: str, keys: list, psnr: list, ssim: list):
    """顯示 PSNR 和 SSIM 的折線圖"""
    plt.clf()
    plt.figure(figsize=(12, 8))

    # host = HostAxes(fig, [0.1, 0.1, 0.8, 0.8])  # 主座標軸
    # fig.add_axes(host)

    # par = ParasiteAxes(host, sharex=host)  # 共享 x 軸
    # fig.add_axes(par)

    host = host_subplot(111)
    par = host.twinx()

    plt.title("PSNR and SSIM vs Block Size")
    host.set_xlabel("Block Size (K)")
    host.set_ylabel("PSNR Value (dB)")
    par.set_ylabel("SSIM Value")

    (p1,) = host.plot(keys, psnr, "o-", label="PSNR")
    (p2,) = par.plot(keys, ssim, "^-", label="SSIM")

    host.legend()

    host.set_xlim(0, 64)
    host.xaxis.set_major_locator(MultipleLocator(8))
    host.xaxis.set_minor_locator(MultipleLocator(2))

    host.set_ylim(0, 40)
    host.yaxis.set_major_locator(MultipleLocator(5))
    host.yaxis.set_minor_locator(MultipleLocator(1.25))

    par.set_ylim(0, 1)
    par.yaxis.set_major_locator(MultipleLocator(0.125))

    host.grid(True)
    # host.grid(True, "minor", linestyle="--")
    # par.grid(True)

    host.tick_params(axis="both", which="both", direction="in", length=6)
    par.tick_params(axis="both", which="both", direction="in", length=6)

    # plt.tight_layout()
    plt.savefig(filename, dpi=300)
    # plt.show()


def show_diff_plot(
    filename: str,
    x: list,
    ylabel: str,
    y1: list,
    label1: str,
    y2: list,
    label2: str,
    y3: list = None,
    label3: str = None,
    xlabel: str = "Block Size (K)",
    *,
    logscale: bool = False,
):
    """顯示 SSIM 差異的折線圖"""
    plt.clf()
    plt.figure(figsize=(12, 8))

    # plt.title("SSIM Difference")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.plot(x, y1, "o-", label=label1)
    plt.plot(x, y2, "^-", label=label2)
    if y3 and label3:
        plt.plot(x, y3, "s-", label=label3)

    plt.xlim(0, 64)
    plt.gca().xaxis.set_major_locator(MultipleLocator(8))
    plt.gca().xaxis.set_minor_locator(MultipleLocator(2))

    # plt.ylim(0, 1)
    if logscale:
        plt.yscale("log")
    else:
        plt.gca().yaxis.set_major_locator(MaxNLocator(8))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))

    plt.grid(True)
    plt.legend()

    plt.savefig(filename, dpi=300)


plt.show()


if __name__ == "__main__":
    x = list(range(1, 101, 1))

    ########################################

    # with open("results_n.json", "r") as f:
    #     results_n = json.load(f)

    # keys = sorted(map(int, results_n.keys()))

    # eeee = []
    # for k, errors in results_n.items():
    #     errors = [abs(float(errors[i - 1])) for i in x]
    #     eeee.append(errors)
    #     # x = list(range(1, len(errors) + 1))
    #     # print(f"n={k}: {errors}")
    #     # plt.plot(x, errors, "o-")
    #     # plt.yscale("log")
    #     # plt.show()

    # show_diff_plot(
    #     "image/h_refinement.png",
    #     x,
    #     "Relative Error",
    #     eeee[0],
    #     "n=2",
    #     eeee[1],
    #     "n=3",
    #     eeee[2],
    #     "n=4",
    #     xlabel="Mesh (M)",
    #     logscale=True,
    # )

    ########################################

    with open("results_m.json", "r") as f:
        results_m = json.load(f)

    keys = sorted(map(int, results_m.keys()))

    eeee = []
    for k, errors in results_m.items():
        errors = [abs(float(errors[i - 1])) for i in x]
        eeee.append(errors)

    show_diff_plot(
        "image/p_refinement.png",
        x,
        "Relative Error",
        eeee[0],
        "m=4",
        eeee[1],
        "m=16",
        eeee[2],
        "m=64",
        xlabel="Nodes (N)",
        logscale=True,
    )
