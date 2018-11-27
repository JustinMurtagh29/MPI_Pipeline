# This R script performs a one-sample Kolmogorov-Smirnov test against a
# discrete null hypothesis. Its whole purpose it to interface between
# MATLAB and the ks.test function of R's "dgof" (discrete goodness-of-
# fit) package.
#
# Two command-line arguments are expected by this script: the
# alternative hypothesis (see "alternative" option of dgof::ks.test
# function) and the path to a HDF5 file containing the null hypothesis
# and observed data.
#
# The output of this script consists of two numbers: the p-value and the
# statistic of the Kolmogorov-Smirnov test. The two numbers are written
# to the standard output and separated by a space.
#
# This function is called from MATLAB by
# +connectEM/+Specificity/kolmogorovSmirnovTest.m.
#
# The dgof package was written by Taylor B. Arnold and John W. Emerson,
# and was published in Arnold, Emerson (2011). Nonparametric goodness-
# of-fit tests for discrete null distributions. The R Journal.
#
# For more information on the dgof package, please consult its CRAN
# page: https://cran.r-project.org/package=dgof
#
# Written by
#   Alessandro Motta <alessandro.motta@brain.mpg.de>
args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 2) {
    stop("Too few input arguments")
}

alternative <- args[1]
file_path <- args[2]

file <- h5::h5file(file_path, 'r')

null_knots <- file["null_knots"][]
null_masses <- file["null_masses"][]
observations <- file["observations"][]

# sanity checks
if(is.unsorted(null_knots)) {
    stop("Knots of null distribution must be sorted in ascending order")
}

if(length(null_knots) != length(null_masses)) {
    stop("Knots and masses of null distribution must have same number of elements")
}

null_ecdf <- cumsum(c(0, null_masses))
null_ecdf <- null_ecdf / null_ecdf[length(null_ecdf)]
null_ecdf <- stats::stepfun(null_knots, null_ecdf)

ks <- dgof::ks.test(observations, null_ecdf, alternative=alternative)
write(c(ks$p.value, ks$statistic), stdout())
