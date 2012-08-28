/* test_parallel.cxx */

#include <mpi.h>
#include <iostream>
#include <douban/matvec/proxy.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/matvec/vec_container.hpp>

int main(int argc, char *argv[]) {

	int rank, nprocs;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int m = 10, n = 4, sub_n = n / 2;
	int load = m / nprocs; 
	int remainder = m % nprocs; 
	int max_load = load + remainder;
	int *rcounts = new int [nprocs];
	int *displs = new int [nprocs];

	douban::mat_container<double> A(m, n);
	douban::vec_container<double> x(sub_n);
	douban::vec_container<double> y(m);

	A = douban::mat_random(A.dim0(), A.dim1());
	x = douban::vec_random(x.size());

	auto sub_A = mat_cols(A, 0, sub_n);	
	douban::mat_container<double> pA(max_load, sub_n);
	douban::vec_container<double> y_tmp(max_load);


	// Copy value of mat A by each procs
	for(int i = 0; i < load; ++i)
		for(int j = 0; j < sub_n; ++j)
			pA.get(i, j) = A.get(i + rank * load, j);	
	if(remainder != 0 && rank == nprocs - 1)
		for(int i = load; i < max_load; ++i)
			for(int j = 0; j < sub_n; ++j)
				pA.get(i, j) = A.get(i + rank * load, j);


	// General matrix vector multiplication 
	y_tmp = douban::gemv(pA, x);
	
	
	// prepare for MPI_Gatherv
	for(int i = 0; i < nprocs; ++i) {
		rcounts[i] = load;
		displs[i] = i * load;
	}
	if(remainder != 0)
		rcounts[nprocs - 1] = max_load;
	if(rank == nprocs - 1)
		load = max_load;
	
	
	// MPI_Gatherv
	MPI_Gatherv(&y_tmp[0], load, MPI_DOUBLE, &y[0], &rcounts[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	
	// check out the result vec
	if(rank == 0)
		for(int i = 0; i < m; ++i)
			std::cout << "result" << y.get(i) << std::endl;


	delete [] rcounts;
	delete [] displs;
	
	MPI_Finalize();

	return 0;
}
