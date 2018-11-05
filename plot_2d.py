
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
        params = parser.parse_args()
        return params

    args = parse_command_line()
    iw = 0 if args.weighted else 1

    # Read in file
    res = lio.load(args.input)
    n_cosmos = res['n_cosmos']
    bins = res['s']
    s = (bins[:-1]+bins[1:])/2
    s = s.reshape(-1, 1)
    norm = res['norm']
    dist = res['2d']

    # get RR, DR, DD, tpcf and tpcfss
    i = 0
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
    ss = s[:, None] * s[None, :]
    ss = ss.reshape(1, ss.shape[0], ss.shape[1])
    tpcfss = tpcf * ss
    tpcfss_err = tpcf_err * ss

    # create figure and subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    ax1, ax2 = axes

    # tpcf
    im1 = axes[0].imshow(tpcf[iw], origin='lower', interpolation='nearest',
                         extent=[0, 200, 0, 200], vmin=-1.0, vmax=1.0,
                         cmap='viridis')
    fig.colorbar(im1, ax=axes[0], label=r'$\xi(s)$', fraction=0.05)

    # tpcf s**2
    im2 = axes[1].imshow(tpcfss[iw], origin='lower', interpolation='nearest',
                         extent=[0, 200, 0, 200], vmin=0., vmax=200.0,
                         cmap='viridis')
    fig.colorbar(im2, ax=axes[1], label=r'$\xi(s)$', fraction=0.05)

    for ax in axes:
        ax.set(xlabel=r'$\sigma$', ylabel=r'$\pi$')
    fig.tight_layout()

    # Save plot or show plot
    if args.output is not None:
        plt.savefig('{}'.format(args.output), bbox_inches='tight')
    else:
        plt.show()

    print('')