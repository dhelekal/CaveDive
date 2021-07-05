#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

List extract_partiton_times(IntegerVector div_idx, 
                            NumericVector div_times,
                            NumericVector vertex_times_ord,
                            DataFrame node_bounds,
                            DataFrame tip_bounds) {
    List sam_times = List::create();
    List coal_times = List::create();
    NumericVector counts (size(div_idx), 0);
    IntegerVector div_from (size(div_idx), NA_INTEGER);

    bool empty_tips=T;

    for (i=0; i<size(div_idx); ++i) {
        int node_hi = node_bounds["hi"](i);
        int node_lo = node_bounds["lo"](i);
        int tip_hi = tip_bounds["hi"](i);
        int tip_lo = tip_bounds["lo"](i);

        auto node_subset_bounds = std::vector<std::pair<int,int>>();
        auto tip_subset_bounds = std::vector<std::pair<int,int>>();

        NumericVector leaf_times = NumericVector::create();

        for (int j=i; j < size(div_idx); ++j) {
            auto tip_sub_bound = std::pair<int,int>(tip_bounds["lo"](j),tip_bounds["hi"](j));
            if (tip_sub_bound.first() > tip_lo &&
                tip_sub_bound.last() < tip.hi &&
                IntegerVector::is_na(div_from(j))) {

                div_from[j] = i;
                leaf_times.push_back(div_times(j));

                node_subset_bounds.push_back(std::pair<int,int>(node_bounds["lo"](j), node_bounds["hi"](j)));
                tip_subset_bounds.push_back(std::pair<int,int>(tip_bounds["lo"](j), tip_bounds["hi"](j)));

                counts(i)--;
            }

            int coal_exl_sum  = 0; 
            int sam_excl_sum = 0;

            for (int j=0; j < size(tip_subset_bounds); ++j) {
                sam_exl_sum += tip_subset_bounds[j].last() - tip_subset_bounds[j].first() + 1;
                coal_exl_sum += node_subset_bounds[j].last() - node_subset_bounds[j].first() + 1;
            }


            NumericVector p_sam(tip_hi-tip_lo+1 - sam_exl_sum,0); 
            NumericVector p_coal(node_hi-node_lo+1 - coal_exl_sum,0);
            int a_idx <- 0;
            int s_idx <- tip_lo;

            if (sam_excl_sum > 0) {
                for(int j=0; j < size(tip_subset_bounds); ++j){
                    int seg_len = tip_subset_bounds[j].first()-s_idx;
                    if(seg_len > 0) {
                        std::memcpy(p_sam.begin()+a_idx, vertex_times_ord.begin()+(tip_subset_bounds[j].first()-1), (seg_len)*sizeof(double));
                        a_idx = tip_subset_bounds[j].last();
                    }
                }
            }  else {
                std::memcpy(p_sam.begin(), vertex_times_ord.begin()+(tip_lo-1), (tip_hi-tip_lo+1)*sizeof(double));
            }
        }
    }


    List out = List::create(Named("sam_times")=sam_times, 
                            Named("coal_times")=coal_times, 
                            Named("counts")=counts,
                            Named("empty_tips")=empty_tips)
}