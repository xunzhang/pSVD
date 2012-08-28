/* test_parallel.cxx */

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

	int m = 10, n = 4, sub_n = n / 2;
	int load = m / nprocs; // 2
	int  remainder = m % nprocs; // 2
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
	//std::cout << "pA.dim0 is" << pA.dim0() << " pA.dim1 is" << pA.dim1() << std::endl;
	douban::vec_container<double> y_tmp(max_load);

	for(int i = 0; i < load; ++i)
		for(int j = 0; j < sub_n; ++j)
			pA.get(i, j) = A.get(i + rank * load, j);	

	if(remainder != 0) {
		if(rank == (nprocs - 1))
			for(int i = load; i < max_load; ++i)
				for(int j = 0; j < sub_n; ++j)
					pA.get(i, j) = A.get(i + rank * load, j);
	} else {
		for(int i = load; i < max_load; ++i)
			for(int j = 0; j < sub_n; ++j)
				pA.get(i, j) = 0.0;
	}

	y_tmp = douban::gemv(pA, x);

	for(int i = 0; i < max_load; ++i)
		std::cout << "y_tmp(i)" << y_tmp.get(i) << std::endl;

	if(rank == (nprocs - 1))
		load = max_load;

	for(int i = 0; i < nprocs; ++i) {
		rcounts[i] = load;
		displs[i] = i * load;
	}
	if(remainder != 0)
		rcounts[nprocs - 1] = max_load;

	//std::cout << "nprocs is" << nprocs << std::endl;
	std::cout << "load is" << load << std::endl;	
	if(rank == 0) {
		for(int i = 0; i < nprocs; ++i) {
			std::cout << "rcounts is" << rcounts[i] << std::endl;
			std::cout << "displs is" << displs[i] << std::endl;
		}
	}
	std::cout << "ssssssssssssssssss int" << sizeof(int) << std::endl;
	std::cout << "ssssssssssssssssss size_t" << sizeof(size_t) << std::endl;
	std::cout << "load is" << load << std::endl;

	//MPI_Gather(&y_tmp[0], (int)load, MPI_DOUBLE, &y[0], (int)load, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(&y_tmp[0], load, MPI_DOUBLE, &y[0], &rcounts[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);


	for(int i = 0; i < m; ++i)
		for(int j = 0; j < sub_n; ++j)
			std::cout << "sub_A(i, j) is " << sub_A.get(i, j) << std::endl;
	/*
	   std::cout << "hello world!" << std::endl;


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
	   */	
	for(int i = 0; i < sub_n; ++i) 
		std::cout << "x(i) is" << x.get(i) << std::endl;

	if(rank == 0)
		for(int i = 0; i < m; ++i)
			std::cout << "result" << y.get(i) << std::endl;


	delete [] displs;
	MPI_Finalize();
	return 0;
}
