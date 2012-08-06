/* test_svd_tr_demo.cxx */
/* Usage: mpic++ -I/home/wuhong_intern/local/include -L/home/wuhong_intern/local/lib -std=c++0x -std=gnu++0x -lmpi -fopenmp -I/usr/include/mysql -Wall -O0 test_svd_tr_demo.cxx -llapack -lgfortran -lblas -lcblas
*/

#include <vector>
#include <iostream>

#include <douban/matvec/vec_kahan_gemv.hpp>
#include <douban/matvec/vec_ge.hpp>
#include <douban/matvec.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/linalg/svd_tr.hpp>
#include <douban/linalg/svd.hpp>
#include <douban/clog.hpp>

int main() 
{
	//double Av_vec[24] = {20, 30, 9, 7, 40, 10, 14, 10, 50, 11, 13, 18, 11, 60, 19, 14, 70, 16, 12, 17, 15, 80, 18, 90}; // size is nnz
	double Av_vec[24] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
	int Ai_vec[24] = {0, 1, 4, 1, 2, 3, 0, 2, 3, 4, 6, 7, 3, 4, 6, 3, 5, 6, 0, 2, 5, 6, 3, 7}; // size is nnz
	int Ap_vec[9] = {0, 1, 3, 6, 12, 15, 18, 22, 24}; // size is m + 1

	int m = 8, n = 8, nnz = 24, k = 5, p = 1;
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

	for(int i = 0; i < 4; i++)
		std::cout << "S(i) is " << S.get(i) << std::endl;

	return 0;
}
