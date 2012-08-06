/* run_svd_tr.cpp */
/* Usage: mpic++ -I/home/wuhong_intern/local/include -L/home/wuhong_intern/local/lib -std=c++0x -std=gnu++0x -lmpi -fopenmp -I/usr/include/mysql -Wall -O0 run_svd_tr.cpp -llapack -lgfortran -lblas -lcblas
*/
#include <mpi.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <douban/clog.hpp>
#include <douban/matvec.hpp>
#include <douban/linalg/svd.hpp>
#include <douban/linalg/svd_tr.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/matvec/vec_kahan_gemv.hpp>

static size_t randint(size_t a, size_t b) {	
	srand((size_t)time(0));
	return size_t(a + (RAND_MAX * rand() + rand()) % (b + 1 - a));
}

int main(int argc, char *argv[]) {
	
	MPI_Init(&argc, &argv);
	
	int res, display = 10;
	double sparsity = 0.01;
	int m = 150000, n = 150000, nnz = 0, k = 200, p = 200;

	if(nnz == 0) 
		nnz = randint((size_t)1, (size_t)m * n * sparsity);
	
	std::cout << "nnz is " << nnz << std::endl;

	if(k > std::min(m, n)) 
		return 1;
	
	std::vector<double> Av(nnz);
	std::vector<size_t> Ap(m + 1), Ai(nnz);
	douban::csr_random(m, n, nnz, Av, Ap, Ai);

	auto A = douban::make_mat_csr(m, n, &Av, &Ap, &Ai);	
	//for(int i = 0; i < m; i++)
	//	for(int j = 0; j < n; j++)
	//		std::cout << "matrix value is " << A.get(i,j) << std::endl;
	douban::mat_container<double> U(m, k);
	douban::mat_container<double> V(n, k);
	douban::vec_container<double> S(k);

	res = douban::linalg::svd_tr(A, U, S, V, p, 1.e-7, display);
	
	std::cout << "U.dim0 is" << U.dim0() << " U.dim1 is" << U.dim1() << std::endl;
	std::cout << "S.size is " << S.size() << std::endl;
	std::cout << "V.dim0 is" << V.dim0() << " V.dim1 is" << V.dim1() << std::endl;
	std::cout << "The result is " << res << std::endl;

	for(int i = 0; i < 5; i++)
		std::cout << "S(i) is " << S.get(i) << std::endl;

	MPI_Finalize();
}
