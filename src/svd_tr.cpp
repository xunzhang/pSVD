/* svd_tr.cpp */

#include "svd_tr.hpp"
#include "mat_csr.hpp"
#include "mat_container.hpp"

#include <vector>
#include <algorithm>

using std::vector;

int svd_tr(mat_csr &A,
	    mat_container<double> &U,
	    vector<double> &S,
	    mat_container<double> &V,
	    size_t over_sample,
	    double eps) {
	
	size_t m = U.dim0(), n = V.dim0();
	size_t nev = std::max(U.dim1(), std::max(V.dim1(), S.size()));
	if(nev <= 0) return 0;
	if(nev + (over_sample ? over_sample : nev) >= std::min(m, n)) {
		return svd_ge(A, U, S, V);
	}
	return svd_tr();
}

