# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#'Compute likelihood of a realisation of a homogeneous coalescent.
#' 
#' @param time_intervals waiting times.
#' @param lineage_count number of lineages for each waiting time.
#' @param pop_size effective population size \code{Neg}.
#' @return the log likelihood for the given parameters.
#' @export
coalescent_loglh <- function(sampling_times, coalescent_times, pop_size, t_max) {
    .Call('_CaveDive_coalescent_loglh', PACKAGE = 'CaveDive', sampling_times, coalescent_times, pop_size, t_max)
}

#'Compute likelihood of a realisation of an exponential growth coalescent.
#' 
#' @param sampling_times times of leaves.
#' @param coalescent_times times of coalescent events.
#' @param lambda growth rate.
#' @param pop_size final effective pop_size at time of most recent sample.
#' @return the log likelihood for the given parameters.
#' @export
exponential_coalescent_loglh <- function(sampling_times, coalescent_times, lambda, pop_size) {
    .Call('_CaveDive_exponential_coalescent_loglh', PACKAGE = 'CaveDive', sampling_times, coalescent_times, lambda, pop_size)
}

#'Compute likelihood of a realisation of an half/logistic growth coalescent.
#' 
#' @param sampling_times times of leaves.
#' @param coalescent_times times of coalescent events.
#' @param t0 time of divergence.
#' @param growth rate.
#' @param N asymptotic population size.
#' @return the log likelihood for the given parameters.
#' @export
logexp_coalescent_loglh <- function(sampling_times, coalescent_times, t0, r, N, t_max) {
    .Call('_CaveDive_logexp_coalescent_loglh', PACKAGE = 'CaveDive', sampling_times, coalescent_times, t0, r, N, t_max)
}

#'Compute likelihood of a realisation of an half/logistic growth coalescent.
#' 
#' @param sampling_times times of leaves.
#' @param coalescent_times times of coalescent events.
#' @param t0 time of divergence.
#' @param growth rate.
#' @param N asymptotic population size.
#' @return the log likelihood for the given parameters.
#' @export
sat_coalescent_loglh <- function(sampling_times, coalescent_times, t0, r, N, t_max) {
    .Call('_CaveDive_sat_coalescent_loglh', PACKAGE = 'CaveDive', sampling_times, coalescent_times, t0, r, N, t_max)
}

extract_partition_times_fast <- function(div_idx, div_times, vertex_times_ord_nodes, vertex_times_ord_tips, node_bounds, tip_bounds) {
    .Call('_CaveDive_extract_partition_times_fast', PACKAGE = 'CaveDive', div_idx, div_times, vertex_times_ord_nodes, vertex_times_ord_tips, node_bounds, tip_bounds)
}

