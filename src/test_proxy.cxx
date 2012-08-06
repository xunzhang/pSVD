/* test_svd_tr_large.cxx */
/* Usage: mpic++ -I/home/wuhong_intern/local/include -L/home/wuhong_intern/local/lib -std=c++0x -std=gnu++0x -lmpi -fopenmp -I/usr/include/mysql -Wall -O0 test_proxy.cxx -llapack -lgfortran -lblas -lcblas
*/

#include <douban/matvec/proxy.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/matvec/vec_container.hpp>


int main(void)
{
	size_t m = 10, n = 8;
	douban::mat_container<double> A(m, n);
	douban::vec_container<double> x(n), y(m);

	A = douban::mat_random(A.dim0(), A.dim1());

	y = douban::vec_random(y.size());

	x = douban::vec_random(x.size());

	douban::vec_container<double> residual(y.size());
	douban::vec_container<double> step(x.size());

	double aa = A.square_sum();

	for(int i = 0; i < 10; i++) {
		residual = douban::gemv(A, x) - y;
		step = (2 / aa) * douban::gemv(A.trans(), residual);
		x += step;
	}

	return 0;
}
