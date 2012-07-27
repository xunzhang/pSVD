/* test_parallel.cxx */
/*
 * mpic++ -I/home/wuhong_intern/local/include -L/home/wuhong_intern/local/lib -std=c++0x -std=gnu++0x -lmpi -fopenmp -I/usr/include/mysql -Wall -O0 test_dense_mv_parallel_s1.cxx -llapack -lgfortran -lblas -lcblas
 *
 *
 * Using douban::mat_container<double> for parallel conputing
 *
 */

#include <iostream>
#include <douban/matvec/proxy.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/matvec/vec_container.hpp>
#include "douban/mpi.hpp"
#include <mpi.h>

int main(int argc, char *argv[]) {
	int rank, nprocs;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	size_t m = 4, n = 2;
	douban::mat_container<double> A(m, n);
	douban::vec_container<double> x(n);
	//if(rank == 0)
	douban::vec_container<double> y(m * nprocs);
	douban::vec_container<double> res(m);

	A = douban::mat_random(A.dim0(), A.dim1());
	x = douban::vec_random(x.size());

	int load = m / nprocs;

	douban::mat_container<double> pA(m, n);
	douban::vec_container<double> y_tmp(m);

	for(int i = 0; i < m; ++i)
		for(int j = 0; j < n; ++j)
			pA.get(i, j) = 0.0;
	for(int i = 0; i < load; ++i)
		for(int j = 0; j < n; ++j)
			pA.get(i + rank * load, j) = A.get(i + rank * load, j);	

	std::cout << "pA.dim0 is" << pA.dim0() << " pA.dim1 is" << pA.dim1() << std::endl;

	y_tmp = douban::gemv(pA, x);

	MPI_Gather(&y_tmp[0], m, MPI_DOUBLE, &y[0], m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank == 0) {
		for(int i = 0; i < m; ++i)
			for(int j = 1; j < nprocs; ++j)
				y.get(i) += y.get(j * m + i);
		for(int i = 0; i < m; ++i)
			res.get(i) = y.get(i);
	}	

	/*
	   std::cout << "hello world!" << std::endl;

	   for(int i = 0; i < load; ++i)
	   std::cout << "y_tmp(i)" << y_tmp.get(i) << std::endl;


	   for(int i = 0; i < m; ++i) 
	   for(int j = 0; j < n; ++j)
	   std::cout << "A(i, j) is" << A.get(i, j) << std::endl;
	   if(rank == 0)	
	   for(int i = 0; i < load; ++i) 
	   for(int j = 0; j < n; ++j)
	   std::cout << "pA(i, j) is" << pA.get(i, j) << std::endl;
	   if(rank == 1)	
	   for(int i = 0; i < load; ++i) 
	   for(int j = 0; j < n; ++j)
	   std::cout << "pA(i, j) is" << pA.get(i, j) << std::endl;

	   for(int i = 0; i < n; ++i) 
	   std::cout << "x(i) is" << x.get(i) << std::endl;
*/
	if(rank == 0) 
		for(int i = 0; i < m; ++i)
			std::cout << "" << y.get(i) << std::endl;
	   

	MPI_Finalize();
	return 0;
}
