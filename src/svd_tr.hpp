/* svd_tr.hpp */
#ifndef SVD_TR_HPP
#define SVD_TR_HPP

#include <vector>
#include "mat_csr.hpp"
#include "mat_container.hpp"

int svd_tr(mat_csr &A, 
            mat_container<double> &U, 
            std::vector<double> &S, 
	    mat_container<double> &V, 
	    size_t over_sample = 0, 
	    double eps = 1.e-7);


#endif
