""" Plotting for NPZ output """

import argparse

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import lib.io as lio


mpl.rc('font', size=15)

if __name__ == '__main__':
    """ Main """

    def parse_command_line():
        parser = argparse.ArgumentParser(description='fast plot')
        parser.add_argument('input',
                            help    = 'input files',
                            type    = str)
        parser.add_argument('--weighted', '-w',
                            help    = 'enable weight',
                            action  = 'store_true',
                            default = False)
        parser.add_argument('--output', '-o',
                            help    = 'output',
                            dest    = 'output',
                            default = None,
                            type    = str)
        parser.add_argument('--scatter', '-s',
                            dest    = 'scatter',
                            action  = 'store_true',
                            default = False)
        parser.add_argument('--error', '-e',
                            action  = 'store_true',
                            default = False)

        params = parser.parse_args()
        return params

    args = parse_command_line()

    # create figure and subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    # Read in file
    tpcf_list = lio.load(args.input)
    iw = 0 if args.weighted else 1

    for i, corr in enumerate(tpcf_list):

        # get RR, DR, DD, tpcf and tpcfss
        rr = corr.rr[iw]
        dr = corr.d1r2[iw]
        dd = corr.dd[iw]
        tpcf = corr.tpcf[iw]
        tpcfss = corr.tpcfss[iw]
        bins = corr.bins

        # Plot
        label = '{}'.format(i+1)
        if not args.error:
            if not args.scatter:
                axes[0, 0].hist(
                    bins[:-1], bins=bins, weights=rr,
                    histtype='step', label=label)
                axes[0, 1].hist(
                    bins[:-1], bins=bins, weights=dr,
                    histtype='step', label=label)
                axes[0, 2].hist(
                    bins[:-1], bins=bins, weights=dd,
                    histtype='step', label=label)
                axes[1, 0].hist(bins[:-1], bins=bins, weights=tpcf,
                                histtype="step", label=label)
                axes[1, 1].hist(bins[:-1], bins=bins, weights=tpcfss,
                                histtype="step", label=label)
            else:
                s = (bins[:-1]+bins[1:])/2
                axes[0, 0].scatter(s, rr, label=label, marker='.')
                axes[0, 1].scatter(s, dr, label=label, marker='.')
                axes[0, 2].scatter(s, dd, label=label, marker='.')
                axes[1, 0].scatter(s, tpcf, label=label, marker='.')
                axes[1, 1].scatter(s, tpcfss, label=label, marker='.')
        else:
            rr_err = corr.rr_err[iw]
            dr_err = corr.d1r2_err[iw]
            dd_err = corr.dd_err[iw]
            tpcf_err = corr.tpcf_err[iw]
            tpcfss_err = corr.tpcfss_err[iw]

            s = (bins[:-1]+bins[1:])/2
            xerr = (s[1]-s[0])/2
            axes[0, 0].errorbar(s, rr, yerr=rr_err, xerr=xerr,
                                label=label, fmt='--.')
            axes[0, 1].errorbar(s, dr, yerr=dr_err, xerr=xerr,
                                label=label, fmt='--.')
            axes[0, 2].errorbar(s, dd, yerr=dd_err, xerr=xerr,
                                label=label, fmt='--.')
            axes[1, 0].errorbar(s, tpcf, yerr=tpcf_err, xerr=xerr,
                                label=label, fmt='--.')
            axes[1, 1].errorbar(s, tpcfss, yerr=tpcfss_err, xerr=xerr,
                                label=label, fmt='--.')
        axes[1, 2].axis('off')

    # set plot labels and legends
    y_label = ["RR(s)", "DR(s)", "DD(s)"] + [r"$\xi(s)$"] + [r"$\xi(s)s^{2}$"]
    for i, ax in enumerate(axes.flat[:-1]):
        ax.set(xlabel="s [Mpc/h]", ylabel=y_label[i], xlim=(0, 200))
        ax.legend()
    fig.tight_layout()

    # Save plot or show plot
    if args.output is not None:
        plt.savefig('{}'.format(args.output), bbox_inches='tight')
    else:
        plt.show()
