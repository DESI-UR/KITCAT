""" Get boundaries of given catalogs """

# Standard Python modules
import sys

# Python modules
import numpy
from astropy.table import Table

# User-defined modules
from lib.cosmology import Cosmology

RAD2DEG = 180./numpy.pi

def main():
    """ Finding the absolute limit of each catalog """
    fname_list = sys.argv[1:-1]
    s_max = float(sys.argv[-1])

    # Sentinel values
    ra_min = 100000
    ra_max = -100000
    dec_min = 100000
    dec_max = -100000
    z_min = 100000
    z_max = -100000

    # Loop over all catalog
    for i, fname in enumerate(fname_list):

       # open .fits file
        table = Table.read(fname)
        dec = table['dec'].data
        ra = table['ra'].data
        z = table['z'].data

        print("Catalog {}".format(i))
        print("RA:  [{}, {}]".format(ra.min(), ra.max()))
        print("DEC: [{}, {}]".format(dec.min(), dec.max()))
        print("Z:   [{}, {}]".format(z.min(), z.max()))

        ra_min = min(ra.min(), ra_min)
        ra_max = max(ra.max(), ra_max)
        dec_min = min(dec.min(), dec_min)
        dec_max = max(dec.max(), dec_max)
        z_min = min(z.min(), z_min)
        z_max = max(z.max(), z_max)

    print("All:")
    print("RA:  [{}, {}]".format(ra_min, ra_max))
    print("DEC: [{}, {}]".format(dec_min, dec_max))
    print("Z:   [{}, {}]".format(z_min, z_max))

if __name__ == "__main__":
    main()
