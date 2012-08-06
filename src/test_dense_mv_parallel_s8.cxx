/* test_parallel.cxx */

#include <mpi.h>
#include <iostream>
#include <douban/matvec/proxy.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/matvec/vec_container.hpp>

template<class T, class P> void sub_task(T& mat, P& vec, P& res) {
	
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	int m = mat.dim0(), n = mat.dim1();
	int load = m / nprocs;
	int remainder = m % nprocs;
	int max_load = load + remainder;
	int *rcounts = new int [nprocs];
	int *displs = new int [nprocs];
	
	douban::vec_container<double> x_tmp(n);
	douban::mat_container<double> pA(max_load, n);
	douban::vec_container<double> y_tmp(max_load);

	
	// copy value of vec by each procs
	for(int i = 0; i < n; ++i)
		x_tmp.get(i) = vec.get(i);


	// copy value of mat by each procs
	for(int i = 0; i < load; ++i)
		for(int j = 0; j < n; ++j)
			pA.get(i, j) = mat.get(i + rank * load, j);
	if(remainder != 0 && rank == nprocs - 1)
		for(int i = load; i < max_load; ++i)
			for(int j = 0; j < n; ++j)
				pA.get(i, j) = mat.get(i + rank * load, j);


	// General matrix vector multiplication
	y_tmp = douban::gemv(pA, x_tmp);
	
	
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
	MPI_Gatherv(&y_tmp[0], load, MPI_DOUBLE, &res[0], &rcounts[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	delete [] rcounts;
	delete [] displs;
}
int main(int argc, char *argv[]) {
	
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int m = 10, n = 4;
	
	douban::mat_container<double> A(m, n);
	douban::vec_container<double> x(m);
	douban::vec_container<double> y(n/2);
	auto Tj = make_vec(&y, 3);
	
	A = douban::mat_random(A.dim0(), A.dim1());
	x = douban::vec_random(x.size());
	
	auto sub_A = mat_cols(A, 0, n / 2);
	auto At = sub_A.trans();	
	At = sub_A.trans();
	auto xx = A.col(3);

	for(int i = 0; i < m; ++i)
		x.get(i) = xx.get(i);	
	// function call
	sub_task(At, x, y);

	for(int i = 0; i < n / 2; ++i)
		xx.get(i) = y.get(i);

	// check out the result vec
	if(rank == 0)
		for(int i = 0; i < n / 2; ++i)
			std::cout << "result" << y.get(i) << std::endl;
		
	MPI_Finalize();

	return 0;
}
