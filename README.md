echoIBM R package
=====

This R package provides utilities for simulating multibeam sonar data (echosounder, fishery sonar, 3D-sonar). Writes simulated data in the TSD format

Version: 1.1
Required R version: 3.5

Installation
=====

``` r
# Install the packages that echoIBM depends on. Note that this updates all the specified packages to the latest (binary) version. To skip installing already installed packages, run install.packages(setdiff(dep.pck, installed.packages()[,"Package"]), repos="http://cran.us.r-project.org") instead:
dep.pck <- c("devtools", "akima", "ccaPP", "data.table", "fBasics", "fields", "fpc", "gdata", "gsl", "pbapply", "SoDA", "XML")
install.packages(dep.pck, repos="http://cran.us.r-project.org", type="binary")

# Install echoIBM and also the packages that echoIBM depends on which are on GitHub (by Holmin):
# On Windows you will need Rtools to complete the installations.
# Check whether you have Rtools by running Sys.getenv('PATH'),
#   and go to https://cran.r-project.org/bin/windows/Rtools/ to install Rtools if not.
# Be sure to check the box "Add rstools to system PATH" when installing Rtools.
# Note that if you need to run R as administrator due to security settings,
#   it is advised to install the pakcages in plain R, and NOT using Rstudio.
# Close Rstudio, open R and run the installation, and reopen Rstudio.

dep.pck.git <- c("arnejohannesholmin/TSD", "arnejohannesholmin/SimradRaw", "arnejohannesholmin/sonR", "arnejohannesholmin/echoIBM")
# If you want to install the lastest development versions, run devtools::install_github(dep.pck.git, ref="develop") instead:
devtools::install_github(dep.pck.git)

```

# For changes log see https://github.com/arnejohannesholmin/echoIBM/blob/master/NEWS

Examples
=====

``` r
# Examples of simulation using the echoIBM package will be out end of 2017.
```

License
=====

The echoIBM package is licensed under the LGPL-3.)

