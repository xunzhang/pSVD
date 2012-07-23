/* svd_tr.hpp */
#ifndef SVD_TR_HPP
#define SVD_TR_HPP

#include <vector>
#include "mat_csr.hpp"

void svd_tr(mat_csr &A, 
            std::vector< vector<double> > &U, 
            std::vector<double> &S, 
	    std::vector< vector<double> > &V, 
	    size_t over_sample = 0, 
	    double eps = 1.e-7);


#endif
