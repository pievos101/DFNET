% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train.R
\name{learn_decisions}
\alias{learn_decisions}
\title{Decision tree learning from modules}
\usage{
learn_decisions(
  raw_modules,
  features,
  target,
  flatten.sep = "$",
  importance = "impurity",
  splitrule = "gini"
)
}
\arguments{
\item{raw_modules}{list of numeric vectors. The raw modules.}

\item{features}{numeric matrix or 3D array. The features to train on.}

\item{target}{numeric vector. The target to train towards.}

\item{flatten.sep}{string. Separator to use when flattening features.}

\item{importance}{variable importance mode.
See \link{ranger:ranger}{ranger::ranger}.}

\item{splitrule}{Splitting rule.
See \link{ranger:ranger}{ranger::ranger}.}
}
\value{
A list of shape (\code{trees}, \code{modules},
\code{modules.weights}), where \code{modules} are the sorted
\code{raw_modules} with individual weights \code{modules.weights}, and
\code{trees} contains one ranger decision tree per module.
}
\description{
\code{learn_decisions} uses \code{ranger} to perform feature selection with
respect to \code{raw_modules}.
}
\keyword{internal}
