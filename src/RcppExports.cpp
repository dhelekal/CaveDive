// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// coalescent_loglh
double coalescent_loglh(NumericVector sampling_times, NumericVector coalescent_times, const double pop_size, const double t_max);
RcppExport SEXP _CaveDive_coalescent_loglh(SEXP sampling_timesSEXP, SEXP coalescent_timesSEXP, SEXP pop_sizeSEXP, SEXP t_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sampling_times(sampling_timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coalescent_times(coalescent_timesSEXP);
    Rcpp::traits::input_parameter< const double >::type pop_size(pop_sizeSEXP);
    Rcpp::traits::input_parameter< const double >::type t_max(t_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(coalescent_loglh(sampling_times, coalescent_times, pop_size, t_max));
    return rcpp_result_gen;
END_RCPP
}
// exponential_coalescent_loglh
double exponential_coalescent_loglh(NumericVector sampling_times, NumericVector coalescent_times, const double lambda, const double pop_size);
RcppExport SEXP _CaveDive_exponential_coalescent_loglh(SEXP sampling_timesSEXP, SEXP coalescent_timesSEXP, SEXP lambdaSEXP, SEXP pop_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sampling_times(sampling_timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coalescent_times(coalescent_timesSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type pop_size(pop_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(exponential_coalescent_loglh(sampling_times, coalescent_times, lambda, pop_size));
    return rcpp_result_gen;
END_RCPP
}
// logexp_coalescent_loglh
double logexp_coalescent_loglh(NumericVector sampling_times, NumericVector coalescent_times, const double t0, const double r, const double N, const double t_max);
RcppExport SEXP _CaveDive_logexp_coalescent_loglh(SEXP sampling_timesSEXP, SEXP coalescent_timesSEXP, SEXP t0SEXP, SEXP rSEXP, SEXP NSEXP, SEXP t_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sampling_times(sampling_timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coalescent_times(coalescent_timesSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double >::type t_max(t_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(logexp_coalescent_loglh(sampling_times, coalescent_times, t0, r, N, t_max));
    return rcpp_result_gen;
END_RCPP
}
// sat_coalescent_loglh
double sat_coalescent_loglh(NumericVector sampling_times, NumericVector coalescent_times, const double t0, const double r, const double N, const double t_max);
RcppExport SEXP _CaveDive_sat_coalescent_loglh(SEXP sampling_timesSEXP, SEXP coalescent_timesSEXP, SEXP t0SEXP, SEXP rSEXP, SEXP NSEXP, SEXP t_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type sampling_times(sampling_timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coalescent_times(coalescent_timesSEXP);
    Rcpp::traits::input_parameter< const double >::type t0(t0SEXP);
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double >::type t_max(t_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(sat_coalescent_loglh(sampling_times, coalescent_times, t0, r, N, t_max));
    return rcpp_result_gen;
END_RCPP
}
// extract_partition_times_fast
List extract_partition_times_fast(IntegerVector div_idx, NumericVector div_times, NumericVector vertex_times_ord_nodes, NumericVector vertex_times_ord_tips, DataFrame node_bounds, DataFrame tip_bounds);
RcppExport SEXP _CaveDive_extract_partition_times_fast(SEXP div_idxSEXP, SEXP div_timesSEXP, SEXP vertex_times_ord_nodesSEXP, SEXP vertex_times_ord_tipsSEXP, SEXP node_boundsSEXP, SEXP tip_boundsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type div_idx(div_idxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type div_times(div_timesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vertex_times_ord_nodes(vertex_times_ord_nodesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vertex_times_ord_tips(vertex_times_ord_tipsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type node_bounds(node_boundsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type tip_bounds(tip_boundsSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_partition_times_fast(div_idx, div_times, vertex_times_ord_nodes, vertex_times_ord_tips, node_bounds, tip_bounds));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CaveDive_coalescent_loglh", (DL_FUNC) &_CaveDive_coalescent_loglh, 4},
    {"_CaveDive_exponential_coalescent_loglh", (DL_FUNC) &_CaveDive_exponential_coalescent_loglh, 4},
    {"_CaveDive_logexp_coalescent_loglh", (DL_FUNC) &_CaveDive_logexp_coalescent_loglh, 6},
    {"_CaveDive_sat_coalescent_loglh", (DL_FUNC) &_CaveDive_sat_coalescent_loglh, 6},
    {"_CaveDive_extract_partition_times_fast", (DL_FUNC) &_CaveDive_extract_partition_times_fast, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_CaveDive(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
