echoIBM R package
=====

This R package provides utilities for simulating multibeam sonar data (echosounder, fishery sonar, 3D-sonar). Writes simulated data in the TSD format

Version: 1.1
Required R version: 3.3.3

Installation
=====

``` r
# Install the packages that echoIBM depends on. Note that this updates all the specified packages to the latest (binary) version. To skip installing already installed packages, run install.packages(setdiff(dep.pck, installed.packages()[,"Package"]), repos="http://cran.us.r-project.org") instead:
dep.pck <- c("akima", "fBasics", "gdata", "gsl", "pbapply")
install.packages(dep.pck, repos="http://cran.us.r-project.org")

# Install echoIBM and also the packages that echoIBM depends on which are on GitHub (by Holmin):
# On Windows you will need Rtools to complete the installations. Check if you have this by running Sys.getenv('PATH'), and go to https://cran.r-project.org/bin/windows/Rtools/ to install Rtools if not.

dep.pck.git <- c("arnejohannesholmin/TSD", "arnejohannesholmin/SimradRaw", "arnejohannesholmin/sonR", "arnejohannesholmin/echoIBM")
devtools::install_github(dep.pck.git)

```

# For changes log see https://github.com/arnejohannesholmin/echoIBM/NEWS

Examples
=====

``` r
# Examples of simulation using the echoIBM package will be out end of 2017.
```

License
=====

The echoIBM package is licensed under the LGPL-3.)

