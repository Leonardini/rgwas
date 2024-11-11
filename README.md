---
editor_options: 
  markdown: 
    wrap: 72
---

# rgwas

ReverseGWAS

Installation instructions \*\*\* this may need to be moved to a separate
file later

1)  If you do not already have devtools installed, please install it via
    
    install.packages("devtools")
    
    Whether you have it installed or not, start R and run the command
    
    library("devtools")

2)  If you do not already have CPLEX Optimisation Studio installed,
    please install it by following the instructions on the IBM website:
    
    <https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-installing>
    
    Note the installation directory as you may need to provide it in the
    next step.

3)  If you do not already have Rcplex installed, add the following line
    to your .Rprofile (this is normally located in your home directory):
    
    Sys.setenv(CPLEX_BIN="Absolute/Path/To/CPLEX/binary")

For instance, with a standard Mac installation this line typically looks
like:

    Sys.setenv(CPLEX_BIN="/Applications/CPLEX_Studio_Community201/cplex/bin/x86-64_osx/cplex")

4)  Once these steps have been successfully carried out, install rgwas via:
    
    devtools::install_github("bichkd/rgwas")

5)  In order to provide the correct parameter settings to CPLEX, add the
    following line to your .Rprofile (this is normally located in your
    home directory):
    
    Sys.setenv(ILOG_CPLEX_PARAMETER_FILE="Absolute/Path/To/rgwas/MyParameters.prm")
    
    Note that the MyParameters.prm file is distributed alongisde the rgwas package.

6)  To check if everything has been successfully installed, run the commands:

    library(rgwas)
    X = rgwas::mainDriver(inputFile = "InputInitial0NA_3.csv")
