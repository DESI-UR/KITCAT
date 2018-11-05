
import argparse

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import lib.io as lio
import lib.correlation as lcorrelation

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
    iw = 0 if args.weighted else 1

    # create figure and subplots
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    # Read in file
    res = lio.load(args.input)
    n_cosmos = res['n_cosmos']
    bins = res['s']
    s = (bins[:-1]+bins[1:])/2
    s = s.reshape(-1, 1)
    norm = res['norm']
    dist = res['1d']

    for i in range(n_cosmos):

        # get RR, DR, DD, tpcf and tpcfss
        rr = dist['rr'][i]
        rr_err = lcorrelation.get_error(dist['rr'][i])
        dd = dist['dd'][i]
        dd_err = lcorrelation.get_error(dist['dd'][i])
        d1r2 = dist['d1r2'][i]
        d1r2_err = lcorrelation.get_error(dist['d1r2'][i])
        d2r1 = dist['d2r1'][i]
        d2r1_err = lcorrelation.get_error(dist['d2r1'][i])
        tpcf, tpcf_err = lcorrelation.tpcf(
            rr          = rr,
            dd          = dd,
            d1r2        = d1r2,
            d2r1        = d2r1,
            norm_rr     = norm['rr'],
            norm_dd     = norm['dd'],
            norm_d1r2   = norm['d1r2'],
            norm_d2r1   = norm['d2r1'])
        tpcfss = tpcf * s**2
        tpcfss_err = tpcf_err * s**2

        # Plot
        label = '{}'.format(i+1)
        if not args.error:
            if not args.scatter:
                axes[0, 0].hist(
                    bins[:-1], bins=bins, weights=rr[iw],
                    histtype='step', label=label)
                axes[0, 1].hist(
                    bins[:-1], bins=bins, weights=dd[iw],
                    histtype='step', label=label)
                axes[1, 0].hist(
                    bins[:-1], bins=bins, weights=d1r2[iw],
                    histtype='step', label=label)
                axes[1, 1].hist(
                    bins[:-1], bins=bins, weights=d2r1[iw],
                    histtype='step', label=label)
                axes[0, 2].hist(bins[:-1], bins=bins, weights=tpcf[iw],
                                histtype="step", label=label)
                axes[1, 2].hist(bins[:-1], bins=bins, weights=tpcfss[iw],
                                histtype="step", label=label)
            else:
                s = (bins[:-1]+bins[1:])/2
                axes[0, 0].scatter(s, rr[iw], label=label, marker='.')
                axes[0, 1].scatter(s, dd[iw], label=label, marker='.')
                axes[1, 0].scatter(s, d1r2[iw], label=label, marker='.')
                axes[1, 1].scatter(s, d2r1[iw], label=label, marker='.')
                axes[0, 2].scatter(s, tpcf[iw], label=label, marker='.')
                axes[1, 2].scatter(s, tpcfss[iw], label=label, marker='.')
        else:
            xerr = (bins[1:]-bins[:-1])/2
            xerr = xerr.reshape(-1, 1)
            axes[0, 0].errorbar(s, rr[iw], yerr=rr_err[iw], xerr=xerr,
                                label=label, fmt='--.')
            axes[0, 1].errorbar(s, dd[iw], yerr=dd_err[iw], xerr=xerr,
                                label=label, fmt='--.')
            axes[1, 0].errorbar(s, d1r2[iw], yerr=d1r2_err[iw], xerr=xerr,
                                label=label, fmt='--.')
            axes[1, 1].errorbar(s, d1r2[iw], yerr=d1r2_err[iw], xerr=xerr,
                                label=label, fmt='--.')
            axes[0, 2].errorbar(s, tpcf[iw], yerr=tpcf_err[iw], xerr=xerr,
                                label=label, fmt='--.')
            axes[1, 2].errorbar(s, tpcfss[iw], yerr=tpcfss_err[iw], xerr=xerr,
                                label=label, fmt='--.')

    # set plot labels and legends
    y_label = ["RR(s)", "DD(s)", r"$\xi(s)$", "D1R2(s)", "D2R1(s)",r"$\xi(s)s^{2}$"]
    for i, ax in enumerate(axes.flatten()):
        ax.set(xlabel="s [Mpc/h]", ylabel=y_label[i], xlim=(0, 200))
        ax.legend()
    fig.tight_layout()

    # Save plot or show plot
    if args.output is not None:
        plt.savefig('{}'.format(args.output), bbox_inches='tight')
    else:
        plt.show()

    print('')