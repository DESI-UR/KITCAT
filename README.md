[![DOI](https://zenodo.org/badge/135309125.svg)](https://zenodo.org/badge/latestdoi/135309125)

# Fast Calculation for Galaxy Two-Point Correlations

A Python implementation of the algorithm described in [A Computationally Efficient Approach for Calculating Galaxy Two-Point Correlations](https://arxiv.org/pdf/1611.09892.pdf).

<!-- toc -->
- [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installing](#installing)
- [Running](#running)
  * [Preprocess](#preprocess)
  * [Combinatorial](#combinatorial)
  * [Integration](#integration)
  * [Plot Output](#plot_output)
- [Configuration File](#configuration-file)
  * [GALAXY Section](#galaxy-section)
  * [RANDOM Section](#random-section)
  * [LIMIT Section](#limit-section)
  * [NBINS Section](#binwidth-section)
  * [COSMOLOGY Section](#cosmology-section)
- [Contributing](#contributing)
- [Versioning](#versioning)
- [Authors](#authors)
- [License](#license)
<!-- tocstop -->


## Getting Started

Installation and running instructions.

### Requirements
Python requirements:
* [Python](https://www.python.org/) (version >3.5.2)
* [AstroPy](http://www.astropy.org) (version >1.2.1)
* [NumPy](http://www.numpy.org) (version >1.13.1)
* [SciPy](https://github.com/scipy/scipy) (version >0.19.1)
* [Scikit-Learn](http://scikit-learn.org/stable/) (version >0.18.1)
* [Matplotlib](https://matplotlib.org/) (Optional) (version >2.0.0)

The sample configurations file rely on DR9-SDSS BOSS survey:
* [SDSS archive](http://www.sdss3.org/dr9/data_access/)


### Installing
To avoid any possible package conflicts, it is strongly encouraged to create a clean `conda` environment. `conda` installation instructions can be found at this [link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/). The new environment can be created as follows:
```
conda create -n KITCATenv python=3
```

Then, KITCAT conda environment can be loaded with:
```
source activate KITCATenv
```

To install KITCAT in the environment:
```
git clone https://github.com/DESI-UR/KITCAT
cd KITCAT
python setup.py install
```

## Running
This implementation has three stages of running: Preprocess, Divide, and Combine.

### Preprocess
Convert galaxy catalog, random catalog, and other parameters (i.e. binnings, cosmological models) into KDTreee, BallTree, and other data structures. These data structures are then stored as binary format into a pickle (.pkl) file, thus further compresses the data. 

If only one cosmological model is specified, then apply the cosmology at the first step and compute comoving distribution instead of redshift distribution.
                     
Options:
    
    - Show help message and exit: 
            -h, --help
    
    - Path to configuration file: 
            -c CONFIG, -C CONFIG, --config CONFIG
    - Output prefix:
            -p PREFIX, -P PREFIX, --prefix PREFIX
    - Total number of Z-slices: 
            -n NSLICE, -N NSLICE, --nslice NSLICE
    - Index of Z-slice. Index runs from 0 to N-1:  
            -i ISLICE, -I ISLICE, --islice ISLICE
    - Set automatic binning: 
            -a, -A, --auto
    - Set binwidth of two-point correlation function. Enable only if auto binning is set: 
            -b BINW, -B BINW, --binwidth BINWIDTH
    - Show program's version number and exit: 
            --version

Example:

The following command: 
```
    KITCAT_preprocess --config=/path/to/sample_conf.cfg --prefix=/path/to/sample_run --islice=0 --nslice=3 
```
will use configuration file at "/path/to/sample_config.cfg". Divide catalogs into 3 redshift slice (z-slice) and process the first slice. The preprocess output will be "sample_run_preprocess.py". 

### Combinatorial
Calculate f(theta), g(theta, z) and DD(s) from preprocess output. This is the most computationally-expensive part of the algorithm. The calculation can be divided into multiple equal-sized child processes. Note: skip the calculation of DD(s) if multiple cosmological models are given in configuration file in PREPROCESS. 

Options:
    
    - Show help message and exit: 
            -h, --help
    - Run prefix:
            -p PREFIX, -P PREFIX, --prefix PREFIX
    - Total number of child processes: 
            -n NJOB, -N NJOB, --nslice NJOB
    - Index of child process. Index runs from 0 to N-1:  
            -i IJOB, -I IJOB, --islice IJOB
    - Save runtime at PREFIX_timer.txt:
            -t, -T, --time
    - Show program's version number and exit: 
            --version

Example:

The following commmand:
```
    KITCAT_combinatorial --prefix=/path/to/sample_run --ijob=0 --njob=10
```
will take preprocess output "/path/to/sample_run_preprocess.pkl". Divide the calculation into 10 child processes and run the first process. The output will be stored at "/path/to/sample_run_divide_000.pkl", with the last three-digit number being the job index).

To perform the calculation without dividing jobs. Simple run
```
    KITCAT_combinatorial --prefix=/path/to/sample_run
```

Following the combinatorial step, user needs to combine the generated files using:
```
    KITCAT_combine --prefix=/path/to/sample_rum
```

### Integration
Perform integration over f(theta), g(theta, r) and P(r) to calculate RR(s), DR(s) and DD(s) (if not already calculated in DIVIDE). Also calculate the two-point correlation function using the Landy-Szalay estimators.

Options:
    
    - Show help message and exit: 
            -h, --help
    - Run prefix:
            -p PREFIX, -P PREFIX, --prefix PREFIX
    - Path to configuration file with cosmological models. If not specified, use the cosmological models from PREPROCESS:
            -c CONFIG, -C CONFIG, --config CONFIG
    - Path to output file with .pkl extension. If not specified, output is saved at PREFIX_output.pkl:
            -o OUTPUT, -O OUTPUT, --output OUTPUT
    - Show program's version number and exit: 
            --version

Example:
The following command
```
    KITCAT_integrate --prefix=/path/to/sample_run --output=/path/to/output 
```
will combine all child processes with prefix "/path/to/sample_run_divide_IJOB.pkl". The output RR(s), DR(s), DD(s), and two-point correlation will be stored at "/path/to/output.pkl".

## Configuration File
This implementation uses Python ConfigParser to read in configuration file. More details on ConfigParser can be found at https://docs.python.org/3/library/configparser.html.

Below are the configuration parameters for each section. Note that configuration parameters are NOT case-sensitive.

### GALAXY and RANDOM Section
Configuration parameters to read galaxy and random catalogs in .fits format.
    
    - PATH: Path to the .fits file that contains the galaxy catalog.
    - DEC: Name of the column that contains declination (DEC).
    - RA: Name of the column that contains right ascension (RA).
    - Z: Name of the column that contains redshift (Z).
    - WEIGHT: Name of the column that contains the total weight of each galaxy.
    - WEIGHT_FKP (optional): Name of the column that contains the FKP weight, a.k.a redshift distribution weight.
    - WEIGHT_SDC (optional): Name of the column that contains the weight due to the distribution of bright stars in the Milky Way.
    - WEIGHT_NOZ (optional): Name of the column that contains the weight due to redshift failure (unable to resolve galaxy redshifts).
    - WEIGHT_CP (optional): Name of the column that contains the weight due to close-pair fiber collision. 

The implementation will search for the key WEIGHT first. If not found, will search for WEIGHT_FKP, WEIGHT_SDC, WEIGHT_NOZ, and WEIGHT_CP and calculate the total weight using the formula: w = w_fkp*w_sdc*(w_noz+w_cp-1)

### NBINS Section
Configuration parameters for the number of bins of the two-point correlation function, and the redshift and angular distribution of the random catalog. Angular variables (RA, DEC, THETA) have unit of degree. Spatial separation have unit of Mpc/h.
    
    - DEC: number of bins for declination for the random angular distribution R(ra, dec).
    - RA: number of bins for right ascension for the random angular distribution R(ra, dec).
    - THETA: binwidth of the pairwise angular distribution f($\theta$) and g($\theta$, r).
    - Z: number of bins for redshift distribution P(z).
    - S: number of bins for the pairwise separation distribution DD(s), DR(s), RR(s), and two-point correlation function.

### LIMIT Section
Configuration parameters for the sample limit of the survey catalog. If the bounds envelope the catalog, then full catalog will be consider. Angular variables (RA, DEC, THETA) have unit of degree. Separation have unit of Mpc/h.
    
    - DEC_MIN: Minimimum declination. 
    - DEC_MAX: Maximimum declination.
    - RA_MIN: Minimum right ascension.
    - RA_MAX: Maximum right ascension.
    - Z_MIN: Minimum redshift.
    - Z_MAX: Maximum redshift.
    - S_MAX: Maximum spatial separation.

### COSMOLOGY Section: 
Cosmological parameters to convert redshift to comoving distance. 
    
    - HUBBLE0: Hubble constant at present (z=0).
    - OMEGA_M0: Relative density of matter at present (z=0).
    - OMEGA_DE0: Relative density of dark energy at present (z=0).

To specify multiple cosmological models, separate the parameters of each model with a comma. For example, the following configuration
    
    - HUBBLE0 = 67, 72, 100
    - OMEGA_M0 = 0.3, 0.4, 0.5
    - OMEGA_DE0 = 0.7 0.6, 0.5

will read in three cosmological models with Hubble constant 67, 72, 100; relative matter density 0.3, 0.4, 0.5; relative dark energy density 0.7, 0.6, 0.5 respectively.

## Contributing
Pending

## Versioning
Pending

## Authors
* **[Tri Nguyen](mailto:tnguy51@u.rochester.edu)**
* **[Tolga Yapici](mailto:tyapici@ur.rochester.edu)**
* **[Regina Demina](http://www.pas.rochester.edu/~regina/)** 
* **[Segev BenZvi](https://www.pas.rochester.edu/~sybenzvi/)**

## License
Pending
