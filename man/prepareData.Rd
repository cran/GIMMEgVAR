% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_prepareGIMMEgVARData.R
\name{prepareData}
\alias{prepareData}
\title{PREPARE DATA FOR FITTING 'GIMMEgVAR'}
\usage{
prepareData(gimmeGVARThreshold = 0.5, RES_matrixData, outputPath)
}
\arguments{
\item{gimmeGVARThreshold}{The cutoff value for group-level paths. Defaults to .50, indicating
that a path must be non-zero across >= .50% of individuals to be
included as a group-level path.}

\item{RES_matrixData}{List containing results of fitting individual 'graphicalVAR' to each individual.}

\item{outputPath}{The user specified path to the directory where results files should be stored.}
}
\description{
Reads in results data from fitting 'graphicalVAR'
and creates the logical matrices for kappa and
beta needed to fit 'GIMMEgVAR'.  The logical matrices determine
which paths will not be regularized in the fitting of 'GIMMEgVAR' to determine
the group level network model.
}
\details{
The following results are returned in the gvarFiles directory:
(1) Returns two data frames that hold the logical matrices for beta and kappa.  These
are named logicalBeta and logicalKappa

(2) Returns two data frames that indicate the proportion of paths across individuals.
These are data frames containing the proportion of paths present across
individuals after fitting 'graphicalVAR'.  These proportions are used to determine
potential group paths which will not be regularized for individuals when fitting the
final 'GIMMEgVAR'. The output files are named proportionKappa.RData and proportionBeta.RData
}
