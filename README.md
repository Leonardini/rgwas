# Installation instructions for ReverseGWAS

## General cross-platform instructions

1)  If you do not already have CPLEX Optimisation Studio installed, please install it by following the instructions on the IBM website:

    <https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-installing>

    Note the installation directory as you may need to provide it in the next step.

2)  If you do not already have Rcplex installed, add the following line to your .Rprofile (this is normally located in your home directory):

```         
    Sys.setenv(CPLEX_BIN="Absolute/Path/To/CPLEX/binary")
```

For instance, with a standard Mac installation this line typically looks like:

```         
    Sys.setenv(CPLEX_BIN="/Applications/CPLEX_Studio_Community201/cplex/bin/x86-64_osx/cplex")
```

3)  If you do not already have devtools installed, please install it from R via

```         
    install.packages("devtools")
```

4)  Once installed (or if you already have it), restart R and run the commands:

```         
    library(devtools)
    devtools::install_github("Leonardini/rgwas")
```

5)  In order to provide the correct parameter settings to CPLEX, add the following line to your .Rprofile (this is normally located in your home directory, and you will need to provide your own path to rgwas):

```         
    Sys.setenv(ILOG_CPLEX_PARAMETER_FILE="Absolute/Path/To/rgwas/MyParameters.prm")
```

Note that the MyParameters.prm file is distributed alongisde the rgwas package.

6)  To check that everything has been successfully installed, run the commands:

```         
    library(rgwas)
    X = rgwas::mainDriver(inputFile = system.file("extdata", "TestInputN5000P5_3.csv", package = "rgwas"), extremeValue = log(1e-3))
    Y = rgwas::optimizingDriver(inputFile = system.file("extdata", "TestInputN5000P5_1.csv", package = "rgwas"), startPValue = log(1e-3))
```

If everything has been configured correctly, all of the following should return TRUE:

```         
    all(X[[1]]$formula == c("", "(p1 OR p3) AND (p1 OR p2)", ""))
    all(Y[[1]]$formula ==       "(p3 OR p5) AND (p1 OR p2)")
    Y[[1]]$LogFisherExactP < X[[1]]$LogFisherExactP[2]
```

7)  If you expect to run this code on large inputs and have patience and a CPLEX license, please try:

```         
    Z = rgwas::validationDriver(inputFile = system.file("extdata", "TestInputN500000P10_3.csv.gz", package = "rgwas"), extremeValue = log(5e-8), shuffle = TRUE)
```

## Specific instructions for using compiled C++ code on Mac or Linux

1)  Begin by following steps 1, 2, 3, and 5 of the cross-platform instructions above, as necessary

2)  If you wish to make use of compiled C++ code, you will need to manually download the rgwas package from GitHub instead of carrying out step 4 of the cross-platform instructions; from the Terminal, run:

```         
    git clone <https://github.com/Leonardini/rgwas.git>
```

3)  On a Mac OS you may need to obtain the GCC compiler by installing homebrew following the instructions at <https://brew.sh/>, then running the following from the Terminal:

```         
    brew install gcc
```

4)  Create a .R subdirectory if one does not already exist in your home directory, and create a Makevars file within that directory if one does not already exist:

```         
    cd ~ 
    mkdir .R 
    cd .R
    touch Makevars
```

5)  Add or append the following two lines to this Makevars file (please note that this step must be done after steps 1 and 2 in the cross-platform instructions, during which the Makevars file must not contain these lines):

```         
    CC = /usr/local/bin/gcc-14
    CXX = /usr/local/bin/g++-14
```

6)  Proceed with a local installation of the rgwas package by changing into the directory where you downloaded rgwas and issuing the command below within R:

```         
    devtools::install_local(path = ".", force = TRUE)
```

7)  To make sure everything worked as expected, please run step 6 (and possibly 7) from the cross-platform instructions above.

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

## Outputs

All exported functions return both the summary statistics as well as the phenotype combination leading to the genotype-phenotype association(s).

Specifically, a typical output includes, for each input genotypic variant, a summary tibble with the following fields:

* SNP: the input genotypic variant             
* time: the time (in seconds) required by the solver, or NA if not recorded           
* gap: the optimality gap at termination (or NA if optimal or not recorded)
* opt: the optimal value of the association statistic found (or NA if infeasible)            
* status: the CPLEX solver status ("optimal", "feasible", "infeasible", "timeout", "nomemory" or "hopeless")          
* formula: the formula representing the phenotype combination       
* TP, FP, FN, TN: the number of true positives, false positives, false negatives, and true negatives           
* numSegments: the number of segments used to represent the $p$-value boundary
* Pos, Neg: the number of positive and negative samples (patients)
* PredPos, PredNeg: the number of positive and negative predictions
* Correct, Incorrect, Total: the number of correct, incorrect, and total predictions
* TPR, TNR: the true positive rate (sensitivity) and true negative rate (specificity)            
* PPV, NPV: the positive and negative predictive values
* Accuracy, OR: the accuracy and odds ratio          
* F1, Jaccard, MCC: the F1 score, Jaccard index, and Matthews correlation coefficient
* LogChiSquaredP, LogFisherExactP: the log of the chi-squared and Fisher exact p-values for the association between the genotype and the combined phenotype
* LogCMHExactP: the log of the Cochrane-Mantel-Haenszel exact p-value, which controls this for the best single phenotype as a binary confounding variable
* bestSingle, bestSingleStat: the single best (most associated) phenotype for the given genotype, and its association statistic 
* computedBounds: if specified, the computed bounds on the association statistic

Optionally, the output may include a phenotype object, which shows the combined phenotype for each patient in the input file.
