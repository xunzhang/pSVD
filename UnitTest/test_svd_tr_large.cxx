/* test_svd_tr_large.cxx */
/* Usage: mpic++ -I/home/wuhong_intern/local/include -L/home/wuhong_intern/local/lib -std=c++0x -std=gnu++0x -lmpi -fopenmp -I/usr/include/mysql -Wall -O0 test_svd_tr_large.cxx -llapack -lgfortran -lblas -lcblas
*/

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <douban/matvec/vec_kahan_gemv.hpp>
//#include <douban/matvec/vec_ge.hpp>
#include <douban/matvec.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/linalg/svd_tr.hpp>
#include <douban/linalg/svd.hpp>
#include <douban/clog.hpp>

static size_t randint(size_t a, size_t b) {	
	srand((size_t)time(0));
	return size_t(a + (RAND_MAX * rand() + rand()) % (b + 1 - a));
}

int main(int argc, char *argv[]) {
	int display = 10;
	double sparsity = 0.01;
	//size_t max_iter = 50;
	//size_t m = 100, n = 100, nnz = 0, k = 50, p = 50;
	int m = 100, n = 100, nnz = 0, k = 10, p = 10;

	if(nnz == 0) 
		nnz = randint((size_t)1, (size_t)m * n * sparsity);
	
	std::cout << "nnz is " << nnz << std::endl;

	if(k > std::min(m, n)) 
		return 1;
	
	std::vector<double> Av(nnz);
	std::vector<size_t> Ap(m + 1), Ai(nnz);
	douban::csr_random(m, n, nnz, Av, Ap, Ai);

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

/*
	auto rU = std::ref(U);
	auto rS = std::ref(S);
	auto rV = std::ref(V);
	douban::linalg::svd_tr(A, rU, rS, rV, p, 1.e-7, display);
*/

	std::cout << "The result is " << res << std::endl;

	for(int i = 0; i < k; i++)
		std::cout << "S(i) is " << S.get(i) << std::endl;

	return 0;

}
