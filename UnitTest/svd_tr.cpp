/* svd_tr.cpp */

#include "svd_tr.hpp"
#include "mat_csr.hpp"
#include "mat_container.hpp"
#include "vec_op.hpp"
#include <vector>
#include <algorithm>

using std::vector;

int svd(mat_coo &A, 
		mat_coo &At,
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


	vector<double> alpha(alpha_buf, alpha_buf + k);
	vector<double> beta(beta_buf, beta_buf + k);
	vector<double> rho(rho_buf, rho_buf + k);
	mat_container<double, false> B(b_buf, k, k);

	mat_container<double, false> P(m, k);
	mat_container<double, false> Q(n, k + 1);
	P = 0.0;
	Q = 0.0;
	vec_random<double>((size_t)(Q.col(0).size()), Q.col(0));
	vec_unit(Q.col(0));
	/*
	 * 	vector<double> test;
	 * 		vec_unit(test);
	 * 			*/
	size_t nconv = 0, l = 0, locked = 0;
	size_t thick_cnt = 0, cnt = 0, mv_cost = 0;
	double ot = 0;
	while(nconv < nev) {
		thick_cnt += 1;
		cnt += 1;
		//bidiag_gkl_restart(locked, l, k, A, At, alpha, beta, rho, P, Q);
		//		//vec_lim_str(alpha, 30);
		//vec_lim_str(veta, 30);		
		mv_cost += k - l;
		B = 0.0;
		//for(size_t i = locked; i < l; ++i) { B(i, i) = s_buf[i]; B(i, l) = rho[i]; }
		//for(size_t i = l; i < k; ++i) { B(i, i) = alpha[i]; }
		//for(size_t i = l; i + 1 < k; ++i) { B(i, i + 1) = beta[i]; }
		mat_container<double, false> BU(bu_buf, k - locked, k - locked);
		mat_container<double, false> BV(bv_buf, k - locked, k - locked);
		vector<double> BS(s_buf + locked, s_buf + k);

		size_t o_nconv = nconv, n_nconv = locked;
		while(n_nconv < k) {
			if(std::abs(rho[n_nconv]) > eps * s_buf[0]) break;
			n_nconv ++;
		}
		if(n_nconv != o_nconv) {
			cnt = 0;
			l = (n_nconv + k) / 2;
			if(l == n_nconv) l = k;
		}
		if(l > locked) {
		}
		Q.col(l) = Q.col(k);
		nconv = n_nconv;
		if(nconv >= nev) break;
		locked = nconv;
	}

	S.assign(s_buf, s_buf + S.size());
	return 0;
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
	}
	mat_coo At = A.trans();
	return svd(A, At, U, S, V, over_sample, eps);
}


