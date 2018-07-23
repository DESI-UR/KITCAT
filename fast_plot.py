""" Plotting for NPZ output """

import argparse

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from lib.myio import load


mpl.rc('font', size=15)

def main():
    """ Main """
    print('Fast plotting modules')

    parser = argparse.ArgumentParser(
        description='Fast plotting for two-point correlation function.')
    parser.add_argument(
        'filenames', type=str, help='Path to two-point correlation function.')
    parser.add_argument('-w', '-W', '--weighted', action='store_true',
                        default=False)
    parser.add_argument('-o', '-O', '--output', type=str, default=None,
                        help='Path to save plot output. Disable showing plot.')
    parser.add_argument('-s', '-S', '--scatter', action='store_true',
                        default=False)
    parser.add_argument('-e', '-E', '--error', action='store_true',
                        default=False)
    parser.add_argument('--version', action='version', version='KITCAT 1.10')
    args = parser.parse_args()

    # Create figure and subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    # Read in file
    correlations = next(load(args.filenames))

    for i, correlation in enumerate(correlations):
        # Get RR, DR, DD, tpcf and tpcfss
        rand_rand = correlation.get_distr('rr', args.weighted)
        data_rand = correlation.get_distr('dr', args.weighted)
        data_data = correlation.get_distr('dd', args.weighted)
        tpcf = correlation.tpcf(args.weighted)
        tpcfss = correlation.tpcfss(args.weighted)
        bins = correlation.bins

        # Plot
        label = '{}'.format(i+1)
        if not args.error:
            if not args.scatter:
                axes[0, 0].hist(
                    bins[:-1], bins=bins, weights=rand_rand[0],
                    histtype='step', label=label)
                axes[0, 1].hist(
                    bins[:-1], bins=bins, weights=data_rand[0],
                    histtype='step', label=label)
                axes[0, 2].hist(
                    bins[:-1], bins=bins, weights=data_data[0],
                    histtype='step', label=label)
                axes[1, 0].hist(bins[:-1], bins=bins, weights=tpcf[0],
                                histtype="step", label=label)
                axes[1, 1].hist(bins[:-1], bins=bins, weights=tpcfss[0],
                                histtype="step", label=label)
            else:
                s = (bins[:-1]+bins[1:])/2
                axes[0, 0].scatter(s, rand_rand[0], label=label, marker='.')
                axes[0, 1].scatter(s, data_rand[0], label=label, marker='.')
                axes[0, 2].scatter(s, data_data[0], label=label, marker='.')
                axes[1, 0].scatter(s, tpcf[0], label=label, marker='.')
                axes[1, 1].scatter(s, tpcfss[0], label=label, marker='.')
        else:
            s = (bins[:-1]+bins[1:])/2
            xerr = (s[1]-s[0])/2
            axes[0, 0].errorbar(s, rand_rand[0], yerr=rand_rand[1], xerr=xerr,
                                label=label, fmt='--.')
            axes[0, 1].errorbar(s, data_rand[0], yerr=data_rand[1], xerr=xerr,
                                label=label, fmt='--.')
            axes[0, 2].errorbar(s, data_data[0], yerr=data_data[1], xerr=xerr,
                                label=label, fmt='--.')
            axes[1, 0].errorbar(s, tpcf[0], yerr=tpcf[1], xerr=xerr,
                                label=label, fmt='--.')
            axes[1, 1].errorbar(s, tpcfss[0], yerr=tpcfss[1], xerr=xerr,
                                label=label, fmt='--.')
        axes[1, 2].axis('off')

    # set plot labels and legends
    y_label = ["RR(s)", "DR(s)", "DD(s)"] + [r"$\xi(s)$"] + [r"$\xi(s)s^{2}$"]
    for i, ax in enumerate(axes.flat[:-1]):
        ax.set(xlabel="s [Mpc/h]", ylabel=y_label[i], xlim=(0, 150))
        if len(args.filenames) > 1:
            if i in np.arange(0, 3):
                ax.legend(loc='upper left')
            else:
                ax.legend()
    axes[1, 1].set(ylim=(-10, 40))

    fig.tight_layout()

    # Save plot or show plot
    if args.output is not None:
        plt.savefig('{}'.format(args.output), bbox_inches='tight')
    else:
        plt.show()


if __name__ == "__main__":
    main()
