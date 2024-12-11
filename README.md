---
editor_options:
  markdown:
    wrap: 72
output: pdf_document
---

# rgwas

ReverseGWAS

Installation instructions

1)  If you do not already have CPLEX Optimisation Studio installed,
    please install it by following the instructions on the IBM website:

    <https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-installing>

    Note the installation directory as you may need to provide it in the
    next step.

2)  If you do not already have Rcplex installed, add the following line
    to your .Rprofile (this is normally located in your home directory):

    Sys.setenv(CPLEX_BIN="Absolute/Path/To/CPLEX/binary")

    For instance, with a standard Mac installation this line typically
    looks like:

    Sys.setenv(CPLEX_BIN="/Applications/CPLEX_Studio_Community201/cplex/bin/x86-64_osx/cplex")

3)  If you do not already have devtools installed, please install it
    from R via

    install.packages("devtools")

    Once installed (or if you already have it), restart R and run the
    commands:

    library(devtools) devtools::install_github("bichkd/rgwas")

4)  In order to provide the correct parameter settings to CPLEX, add the
    following line to your .Rprofile (this is normally located in your
    home directory, and you will need to provide your own path to
    rgwas):

    Sys.setenv(ILOG_CPLEX_PARAMETER_FILE="Absolute/Path/To/rgwas/MyParameters.prm")

    Note that the MyParameters.prm file is distributed alongisde the
    rgwas package.

5)  If you are in a Mac or Linux environment and can make use of
    compiled C++ code, you will need to download the rgwas package from
    GitHub; from the Terminal, run:

    git clone <https://github.com/bichkd/rgwas.git>

6)  To obtain the needed GCC compiler, you may need to install homebrew
    following the instructions at <https://brew.sh/>, then running the
    following from the Terminal:

    brew install gcc

7)  Create a .R subdirectory if one is not already there in your user
    directory, and create a Makevars file within that directory if one
    does not already exist:

    cd \~ ; mkdir .R ; cd .R ; touch Makevars

    Add or append the following two lines to this Makevars file:

    CC = /usr/local/bin/gcc-14 
    
    CXX = /usr/local/bin/g++-14
    
8)  Proceed with a local installation of the rgwas package by changing into the 
    directory where you downloaded rgwas and issuing the command below within R:

    devtools::install_local(path = ".", force = TRUE,
    configure.vars = "CC=/usr/local/bin/gcc-14 CXX=/usr/local/bin/g++-14")

9)  To check that everything has been successfully installed, run the
    commands:

    library(rgwas)

    X = rgwas::mainDriver(inputFile = "TestInputN5000P5_3.csv",
    extremeValue = log(1e-3))

10)  If you expect to run this code on large inputs and have patience and
    a CPLEX license, please try:

    Y = rgwas::validationDriver(inputFile = "TestInputN500000P10_3.csv",
    extremeValue = log(5e-8))
