#include <vector>
#include <iostream>

#include <douban/matvec/vec_kahan_gemv.hpp>
#include <douban/matvec.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/linalg/svd_tr.hpp>
#include <douban/linalg/svd.hpp>
#include <douban/clog.hpp>

int main() {

	double Av_vec[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; // size is nnz
	int Ai_vec[10] = {0, 3, 1, 2, 4, 0, 1, 2, 3, 4}; // size is nnz
	int Ap_vec[7] = {0, 2, 3, 5, 7, 9, 10}; // size is m + 1

	int m = 6, n = 5, nnz = 10, k = 5, p = 2;
	int display = 10;
	// p is m+1!!	
	std::vector<double> Av(Av_vec, Av_vec + nnz);
	std::vector<int> Ap(Ap_vec, Ap_vec + m + 1), Ai(Ai_vec, Ai_vec + nnz);

	//douban::csr_to_canonical(m, n, Av, Ap, Ai, [&](double v){return v > 0;});
	auto A = douban::make_mat_csr(m, n, &Av, &Ap, &Ai);	
	
	//for(int i = 0; i < m; i++)
	//	for(int j = 0; j < n; j++)
	//		std::cout << "matrix value is " << A.get(i,j) << std::endl;

	douban::mat_container<double> U(m, k);
	douban::mat_container<double> V(n, k);
	douban::vec_container<double> S(k);
	
	//for(int i = 0; i < 5; i++) 
	//	std::cout << "S(i) is " << S.get(i) << std::endl;
	
	int res;
	res = douban::linalg::svd_tr(A, U, S, V, p, 1.e-7, display);

	std::cout << "U.dim0 is" << U.dim0() << " U.dim1 is" << U.dim1() << std::endl;
	std::cout << "S.size is " << S.size() << std::endl;
	std::cout << "V.dim0 is" << V.dim0() << " V.dim1 is" << V.dim1() << std::endl;

	std::cout << "The result is " << res << std::endl;

	for(int i = 0; i < k; i++)
		std::cout << "S(i) is " << S.get(i) << std::endl;
	

	return 0;}
