""" Script for combining job results and calculate DD(s), DR(s), and RR(s) """

# Standard Python module
import argparse
import glob

# User-defined module
import lib.io as lio
import lib.correlation as lcorrelation

if __name__ == "__main__":
    """ integration """

    print('')

    def parse_command_line():
        parser = argparse.ArgumentParser(description='preprocess data')
        parser.add_argument('--prefix', '-p',
                            help    = 'output prefix.',
                            dest    = 'prefix',
                            type    = str)
        parser.add_argument('--output', '-o',
                            help    = 'output name',
                            default = None,
                            dest    = 'output',
                            type    = str)
        params = parser.parse_args()
        return params

    args = parse_command_line()

    # read in helper
    print('reading file')
    fname_list = sorted(glob.glob("%s_divide_*.pkl" % args.prefix))
    helper = None
    for i, fname in enumerate(fname_list):
        print("- %s" % fname)
        if i != 0:
            helper.add(lio.load(fname))
            continue
        helper = lio.load(fname)

    # calculate dd, dr, rr
    print('')
    print('calculate DD(s)')
    dd = helper.get_dd()

    print('')
    print('calculate DR(s)')
    dr = helper.get_dr()

    print('')
    print('calculate RR(s)')
    rr = helper.get_rr()

    # calculate tpcf
    print('')
    print('calculate tpcf')
    tpcf_list = []
    for i in range(dd.shape[0]):
        tpcf = lcorrelation.Correlation(
            dd      = dd[i],
            dr      = dr[i],
            rr      = rr[i],
            bins    = helper.bins.bins('s'),
            norm_rr = helper.norm_rr,
            norm_dr = helper.norm_dr,
            norm_dd = helper.norm_dd)
        tpcf_list.append(tpcf)

    # save output
    if args.output is None:
        lio.save('%s_output.pkl' % args.prefix, tpcf_list)
    else:
        if not args.output.endswith('.pkl'):
            args.output += '.pkl'
        lio.save(args.output, tpcf_list)
