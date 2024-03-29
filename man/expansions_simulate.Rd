% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expansions_simulate.R
\name{expansions_simulate}
\alias{expansions_simulate}
\title{Simulate parameters and event times for a genealogy with clonal expansionss}
\usage{
expansions_simulate(priors, sampling_times, concentration, given = list())
}
\arguments{
\item{priors}{List of priors to simulate parameters from. See `standard priors` for details}

\item{sampling_times}{Vector sampling times, entries must be negative or 0}

\item{concentration}{Scalar concentration hyperparameter for dirichlet expansion model}

\item{given}{A list of variables with values given. Supported names: 'n_exp' - number of expansions, 'tip_colours' - tip expansion assignment, N' - Background population size, 'K' - Expansion carrying capacities, 't_mid' - Expansion times to midpoints, 'div_times' - Expansion divergence times, 'div_from' - parent lineages}
}
\value{
list of: `co` - realisation of the prpcess, `params` - the simulated parameters for the process, `coal_log_lh` - the process likelihood, `param_log_lh` - the prior likelihood
}
\description{
Simulate parameters and event times for a genealogy with clonal expansionss
}
