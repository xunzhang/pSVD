/* svd_tr.cpp */

#include "svd_tr.hpp"
#include "mat_csr.hpp"
#include "mat_container.hpp"

#include <vector>
#include <algorithm>

using std::vector;

int svd(mat_coo &A, 
		mat_container<double> &U,
		vector<double> &S,
		mat_container<double> &V,
		size_t over_sample,
		double eps) {

	size_t m = U.dim0(), n = V.dim0();
	size_t nev = std::max(U.dim1(), std::max(V.dim1(), S.size()));
	
	if(nev <= 0) return 0;
	if(nev > std::min(m, n)) return -1;
	if(over_sample == 0) over_sample = nev;
	
	size_t k = nev + over_sample;
	k = std::min(k, std::min(m, n));
	
	std::vector<double> work(std::max(m, n) * k + 3 * k * k + 4 * k, 0);
	double *mat_buf = &work[0];
	double *b_buf = mat_buf + std::max(m, n) * k;
	double *bu_buf = b_buf + k * k;
	double *bv_buf = bu_buf + k * k;
	double *alpha_buf = bv_buf + k * k;
	double *beta_buf = alpha_buf + k;
	double *rho_buf = beta_buf + k;
	double *s_buf = rho_buf + k;
	
	auto alpha = make_vec(alpha_buf, k);
	auto beta = make_vec(beta_buf, k);
	auto rho = make_vec(rho_buf, k);
	auto B = make_mat_col_major(b_buf, k, k);

	mat_container<double> P(m, k);
	mat_container<double> Q(n, k + 1);
	P = 0;
	Q = 0;
	Q.col(0) = vec_random(Q.col(0).size()) - .5;
	vec_unit(Q.col(0));

	size_t nconv = 0, l = 0, locked = 0;
	size_t thick_cnt = 0, cnt = 0, mv_cost = 0;
	double ot = 0;

}

int svd_tr(mat_coo &A,
		mat_container<double> &U,
		vector<double> &S,
		mat_container<double> &V,
		size_t over_sample,
		double eps) {

	size_t m = U.dim0(), n = V.dim0();
	size_t nev = std::max(U.dim1(), std::max(V.dim1(), S.size()));
	if(nev <= 0) return 0;
	if(nev + (over_sample ? over_sample : nev) >= std::min(m, n)) {
		//return svd_ge(A, U, S, V);
	}
	return svd(A, U, S, V, over_sample, eps);
}

