# Installation instructions for ReverseGWAS

## General cross-platform instructions

1)  If you do not already have CPLEX Optimisation Studio installed, please install it by following the instructions on the IBM website:

    <https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-installing>

    Note the installation directory as you may need to provide it in the next step.

2)  If you do not already have Rcplex installed, add the following line to your .Rprofile (this is normally located in your home directory):

    ```R
    Sys.setenv(CPLEX_BIN="Absolute/Path/To/CPLEX/binary")
    ```

    For instance, with a standard Mac installation this line typically looks like:

    ```R
    Sys.setenv(CPLEX_BIN="/Applications/CPLEX_Studio_Community201/cplex/bin/x86-64_osx/cplex")
    ```

    Then install Rcplex as you normally would. Generally, this would be something like:

    ```R
    install.packages("Rcplex")
    ```

    Please make sure Rcplex is successfully installed before proceeding.

3)  If you do not already have devtools installed, please install it from R via

    ```R         
    install.packages("devtools")
    ```

4)  Once installed (or if you already have it), restart R and run the commands:

    ```  R       
    library(devtools)
    devtools::install_github("Leonardini/rgwas")
    ```
    
    This will install the rgwas package.
    For Linux and Mac, we provide a fast C++ implementation for one part of the approach ("tunnel search"). This requires a compiler that supports the quadmath library.
    If this library is not available, rgwas will fall back on a pure R implementation. See below for special instructions on a Mac.

5)  In order to provide the correct parameter settings to CPLEX, add the following line to your .Rprofile (this is normally located in your home directory, and you will need to provide your own path to rgwas):

    ```R         
    Sys.setenv(ILOG_CPLEX_PARAMETER_FILE="Absolute/Path/To/rgwas/MyParameters.prm")
    ```
    
    Note that the MyParameters.prm file is distributed alongisde the rgwas package.

6)  To check that everything has been successfully installed, run the commands:

    ```R         
    library(rgwas)
    X = rgwas::mainDriver(inputFile = system.file("extdata", "TestInputN5000P5_3.csv", package = "rgwas"), extremeValue = log(1e-3))
    Y = rgwas::optimizingDriver(inputFile = system.file("extdata", "TestInputN5000P5_1.csv", package = "rgwas"), startPValue = log(1e-3))
    ```
    If everything has been configured correctly, all of the following should return TRUE:

    ```R         
    all(X[[1]]$formula == c("", "(p1 OR p3) AND (p1 OR p2)", ""))
    all(Y[[1]]$formula ==       "(p3 OR p5) AND (p1 OR p2)")
    Y[[1]]$LogFisherExactP < X[[1]]$LogFisherExactP[2]
    ```

7)  If you expect to run this code on large inputs and have patience and a CPLEX license, please try:

    ```R         
    Z = rgwas::validationDriver(inputFile = system.file("extdata", "TestInputN500000P10_3.csv.gz", package = "rgwas"), extremeValue = log(5e-8), shuffle = TRUE)
    ```

## Specific instructions for using compiled C++ code on Mac

1)  Begin by following steps 1, 2 and 3 of the cross-platform instructions above, as necessary

2)  On a Mac OS you may need to obtain the GCC compiler (which comes with support for the quadmath library) by installing homebrew following the instructions at <https://brew.sh/>, then running the following from the Terminal:

```         
    brew install gcc
```

3)  Create a .R subdirectory if one does not already exist in your home directory, and create a Makevars file within that directory if one does not already exist:

```         
    cd ~ 
    mkdir .R 
    cd .R
    touch Makevars
```

4)  Add or append the following lines to this Makevars file (please note that this step must be done after steps 1 and 2 in the cross-platform instructions, during which the Makevars file must not contain these lines):

```         
    CC = /usr/local/bin/gcc-14
    CXX = /usr/local/bin/g++-14
    CXX11 = /usr/local/bin/g++-14
    CXX14 = /usr/local/bin/g++-14
    CXX17 = /usr/local/bin/g++-14
    CXX20 = /usr/local/bin/g++-14
```

5)  Proceed with steps 4 through 7 of the cross-platform instructions

# Usage instructions for ReverseGWAS

## Input file formats

Please note that all input files must follow the same format: a patient identifier column named ID (typically a character), followed by one or more genotype columns, followed by one or more phenotype columns.

The genotype and phenotype columns can be either numeric or logical, but if numeric, all values must be 0 or 1; note that missing genotype values are allowed, but missing phenotype values are not.

In addition, the filename must specify the number of genotype columns immediately before the extension, which must be .csv[.gz] or .tsv[.gz] - the compression is optional.

Please see the example files provided with the package.

## Exported functions

There are three exported functions:

`mainDriver`, which takes a multi-genotype and multi-phenotype input with a single parameter setting for K (clause number) and L (clause size), and checks the feasibility of a target value of the association statistic.

`optimizingDriver`, which takes a single genotype and multi-phenotype input with a single K and L, and identifies the optimal value of the association statistic.

`validationDriver`, which performs a discovery-validation split of a multi-genotype and multi-phenotype input file and an exploration of multiple K and L values, and identifies the optimal value of the association statistic for each combination of values.

In every case, both the summary statistics as well as the phenotype combination leading to the genotype-phenotype association(s) are returned.
