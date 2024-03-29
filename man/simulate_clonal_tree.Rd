% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expansions_simulate.R
\name{simulate_clonal_tree}
\alias{simulate_clonal_tree}
\title{Simulate an instance of coalescent process with local population structure}
\usage{
simulate_clonal_tree(
  n_exp,
  N,
  K,
  A,
  sampling_times,
  tip_colours,
  div_times,
  div_cols,
  div_from = NA
)
}
\arguments{
\item{n_exp}{number of expansions}

\item{N}{background population size}

\item{K}{carrying capacities}

\item{A}{growth rates}

\item{sampling_times}{sampling times}

\item{tip_colours}{tip clade assignment}

\item{div_times}{clade divergence times}

\item{div_cols}{vector of colours ordered in sequence}

\item{div_from}{(optional) parent populations for individual clades. If NA will be randomised equiprobably.}
}
\value{
list of: `co` - list containing event times, event colour assignment and parents of clonal expansions; log_lh: the log likelihood
}
\description{
Simulate an instance of coalescent process with local population structure
}
