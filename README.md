# BAPS

Bayesian Analysis of Population Structure. For more information, see reference below.

Corander, J. and Marttinen, P. (2006), Bayesian identification of admixture events using multilocus molecular markers. Molecular Ecology, 15: 2833-2843. doi:10.1111/j.1365-294X.2006.02994.x

# Installation

## Running MATLAB from the command line

Open a terminal. From the root folder of this repository, open MATLAB. Assuming MATLAB is installed in `/path/to/matlab`, this can be done by typing:

```bash
/path/to/matlab
```

If you don't know where MATLAB is installed, you can find out by opening the MATLAB GUI and typing:

```MATLAB
matlabroot
```

On Linux and Windows, respectively, this is `/usr/local/MATLAB/R2016b` and `C:\Program Files\MATLAB\R2016b` by default.

Finally, run MATLAB from the command line by typing:

```bash
/path/to/matlab -nodesktop -nosplash
```

## Compiling BAPS

Once MATLAB is running, compile BAPS by typing:

```MATLAB
run add_BAPS_to_path.m
run compileBaps6.m
```

This should create a `BAPS_package` folder in the root directory of this repository.
