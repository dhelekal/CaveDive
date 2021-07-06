#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List extract_partition_times_fast(IntegerVector div_idx, 
    NumericVector div_times,
    NumericVector vertex_times_ord_nodes,
    NumericVector vertex_times_ord_tips,
    DataFrame node_bounds,
    DataFrame tip_bounds) 
{
    List sam_times = List::create();
    List coal_times = List::create();
    IntegerVector counts (div_idx.size(), 0);
    IntegerVector div_from (div_idx.size(), NA_INTEGER);

    bool empty_tips=false;

    IntegerVector node_bounds_col_hi = node_bounds["hi"];
    IntegerVector node_bounds_col_lo = node_bounds["lo"];

    IntegerVector tip_bounds_col_hi = tip_bounds["hi"];
    IntegerVector tip_bounds_col_lo = tip_bounds["lo"];

    for (int i=0; i<div_idx.size(); ++i) 
    {
        int c_idx = div_idx(i)-1;
        int node_hi = node_bounds_col_hi(c_idx);
        int node_lo = node_bounds_col_lo(c_idx);
        int tip_hi = tip_bounds_col_hi(c_idx);
        int tip_lo = tip_bounds_col_lo(c_idx);

        auto node_subset_bounds = std::vector<std::pair<int,int>>();
        auto tip_subset_bounds = std::vector<std::pair<int,int>>();

        NumericVector leaf_times = NumericVector::create();

        for (int j=0; j < i; ++j) 
        {
            int d_idx = div_idx(j)-1;
            auto tip_sub_bound = std::pair<int,int>(tip_bounds_col_lo(d_idx), tip_bounds_col_hi(d_idx));

            if (tip_sub_bound.first >= tip_lo &&
                tip_sub_bound.second <= tip_hi &&
                IntegerVector::is_na(div_from(j))) 
            {   
                div_from[j] = i;
                leaf_times.push_back(div_times(j));

                node_subset_bounds.push_back(std::pair<int,int>(node_bounds_col_lo(d_idx), node_bounds_col_hi(d_idx)));
                tip_subset_bounds.push_back(tip_sub_bound);
            }
        }

        int coal_excl_sum  = 0; 
        int sam_excl_sum = 0;

        for (int j=0; j < tip_subset_bounds.size(); ++j) 
        {
            sam_excl_sum += tip_subset_bounds[j].second - tip_subset_bounds[j].first + 1;
            coal_excl_sum += node_subset_bounds[j].second - node_subset_bounds[j].first + 1;
        }

        std::sort(tip_subset_bounds.begin(), tip_subset_bounds.end(), [](std::pair<int,int> a, std::pair<int,int> b) 
        {
            return a.first < b.first;
        });

        std::sort(node_subset_bounds.begin(), node_subset_bounds.end(), [](std::pair<int,int> a, std::pair<int,int> b) 
        {
            return a.first < b.first;
        });

        counts[i] = tip_hi-tip_lo+1 - sam_excl_sum;
        NumericVector p_sam(counts[i],0.0); 
        NumericVector p_coal(node_hi-node_lo+1 - coal_excl_sum,0.0);
        int a_idx = 0;
        int s_idx = tip_lo;

        if (sam_excl_sum > 0) 
        {
            for(int j=0; j < tip_subset_bounds.size(); ++j)
            {
                int seg_len = tip_subset_bounds[j].first-s_idx;
                if(seg_len > 0) 
                {
                    std::memcpy(p_sam.begin()+a_idx, vertex_times_ord_tips.begin()+(s_idx-1), (seg_len)*sizeof(double));
                    a_idx += seg_len;
                }
                s_idx = tip_subset_bounds[j].second+1;
            }
            int seg_len = tip_hi-s_idx+1;
            if (seg_len > 0) 
            {
                std::memcpy(p_sam.begin()+a_idx, vertex_times_ord_tips.begin()+(s_idx-1), (seg_len)*sizeof(double));
            }
        }  
        else 
        {
            std::memcpy(p_sam.begin()+a_idx, vertex_times_ord_tips.begin()+(tip_lo-1), (tip_hi-tip_lo+1)*sizeof(double));
        }

        a_idx = 0;
        s_idx = node_lo;

        if (sam_excl_sum > 0) 
        {
            for(int j=0; j < node_subset_bounds.size(); ++j)
            {
                int seg_len = node_subset_bounds[j].first-s_idx;
                if(seg_len > 0) 
                {
                    std::memcpy(p_coal.begin()+a_idx, vertex_times_ord_nodes.begin()+(s_idx-1), (seg_len)*sizeof(double));
                    a_idx += seg_len;
                }
                s_idx = node_subset_bounds[j].second+1;
            }
            int seg_len = node_hi-s_idx+1;
            if (seg_len > 0) 
            {
                std::memcpy(p_coal.begin()+a_idx, vertex_times_ord_nodes.begin()+(s_idx-1), (seg_len)*sizeof(double));
            }
        }  
        else 
        {
            std::memcpy(p_coal.begin()+a_idx, vertex_times_ord_nodes.begin()+(node_lo-1), (node_hi-node_lo+1)*sizeof(double));
        }
        if(p_sam.size()==0) 
        {
            empty_tips = true;
        }
        
        for (int j=0; j<leaf_times.size(); ++j) 
        {
            p_sam.push_back(leaf_times[j]);
        }
        std::sort(p_sam.begin(), p_sam.end(), std::greater<double>());
        std::sort(p_coal.begin(), p_coal.end(), std::greater<double>());

        sam_times.push_back(p_sam);
        coal_times.push_back(p_coal);
    }

    List out = List::create(Named("sam.times")=sam_times, 
        Named("coal.times")=coal_times, 
        Named("partition_counts")=counts,
        Named("empty_tips")=empty_tips);
    return(out);
}