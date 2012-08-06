/* run_svd_tr.cpp */

/* Usage: mpic++ -I/home/wuhong_intern/local/include -L/home/wuhong_intern/local/lib -std=c++0x -std=gnu++0x -lmpi -fopenmp -I/usr/include/mysql      -Wall -O0 main.cpp -llapack -lgfortran -lblas -lcblas
*/

#include <map>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <memory.h>
#include <douban/matvec.hpp>
#include <douban/linalg/svd.hpp>
#include <douban/linalg/svd_tr.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/matvec/vec_kahan_gemv.hpp>

int main(int argc, char *argv[]) {

	int rank, nprocs;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int m, n, m_At, n_At, nnz, k = 200, p = 200, display = 10;
	int len, begin, end, s_indx, t_s_indx;
	std::map<int, int> word_kind;
	FILE *fp;


	fp = fopen("pwtk.bin", "rb");
	fread(&m, sizeof(int), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	fread(&nnz, sizeof(int), 1, fp);
	
	/* for demo test
	m = 4;
	n = 3;
	nnz = 6;
	for demo test */

	n_At = m;

	douban::mat_container<double> U(m, k);
	douban::mat_container<double> V(n, k);
	douban::vec_container<double> S(k);

	std::vector<int> ar, ac, Ap, Ai, atr, atc, Atp, Ati;
	std::vector<double> av, Av, Atv;

	/*
	if(rank == 0) {
		std::cout << std::endl;
		std::cout << "Matrix file general info:" << std::endl;
		std::cout << "Dimension size:" << m << " " << n << std::endl;
		std::cout << "Number of nonzeros:" << nnz << std::endl;
	}
	*/
	len = nnz / nprocs;	
	begin = rank * len;
	end = (rank + 1) * len;
	if(rank == nprocs - 1) {
		end = nnz;
		len = nnz - begin;
	}
	std::cout << len << std::endl;

	ar.resize(len);
	ac.resize(len);
	av.resize(len);
	atr.resize(len);
	atc.resize(len);


	fseek(fp, (begin * (sizeof(int)*2 + sizeof(double))), SEEK_CUR);
	for(int i = 0; i < len; ++i) {
		fread(&(ar[i]), sizeof(int), 1, fp);
		fread(&(ac[i]), sizeof(int), 1, fp);
		fread(&(av[i]), sizeof(double), 1, fp);
		atr[i] = ac[i];
		atc[i] = ar[i];
	}
	if(rank == 1) {
		std::cout << "ss" << len << std::endl;
		for(int i = 0; i < 1400000; ++i)
			std::cout << i << std::endl; //<< "nerves "<< atc[i] << std::endl;
	}

	/* for demo test 
	   if(rank == 0) {
	   atc[0] = ar[0] = 0;
	   atc[1] = ar[1] = 0;
	   atc[2] = ar[2] = 1;

	   atr[0] = ac[0] = 0;
	   atr[1] = ac[1] = 1;
	   atr[2] = ac[2] = 0;

	   av[0] = 1;
	   av[1] = 2;
	   av[2] = 3;
	   }
	   if(rank == 1) {
	   atc[0] = ar[0] = 1;
	   atc[1] = ar[1] = 2;
	   atc[2] = ar[2] = 3;

	   atr[0] = ac[0] = 1;
	   atr[1] = ac[1] = 2;
	   atr[2] = ac[2] = 2;

	   av[0] = 4;
	   av[1] = 5;
	   av[2] = 6;
	   }
	   for demo test */

	m = 0;
	for(int i = 0; i < len; ++i) 
		if(ar[i] != ar[i+1])
			m++;
	m_At = 0;
	for(int i = 0; i < len; ++i) {
		if(word_kind.count(atr[i]) == 0) {
			word_kind.insert(std::map<int, int>::value_type(atr[i],1));
			m_At++; 
		}
	} // m_At equals to word_kind.size()

	//std::cout << "m_At" << m_At << std::endl;
	//std::cout << "n_At" << n_At << std::endl;
	//std::cout << "len" << len << std::endl;
	s_indx = ar[0];
	t_s_indx = atr[0];
	//std::cout << "fuck here is" << s_indx << std::endl;
	//std::cout << "fuck here is" << t_s_indx << std::endl;
	if(rank != 0) {
		for(int i = 0; i < len; ++i) ar[i] -= s_indx;
		for(int i = 0; i < len; ++i) atr[i] -= t_s_indx;	
	}
	//std::cout << "m is" << m << std::endl;
	//std::cout << "n is" << n << std::endl;

	Ap.resize(m + 1);
	Av.resize(len);
	Ai.resize(len);

	Atp.resize(m_At + 1);
	Atv.resize(len);
	Ati.resize(len);

	/* for demo test 
	   for(int i = 0; i < 3; ++i) {
	   std::cout << "atr[i] is " << atr[i] << std::endl;
	   std::cout << "atc[i] is " << atc[i] << std::endl;
	   std::cout << "av[i] is " << av[i] << std::endl;
	   }
	   for demo test */

	douban::coo_to_csr(m, n, len, av, ar, ac, Av, Ap, Ai);	
	douban::coo_to_csr(m_At, n_At, len, av, atr, atc, Atv, Atp, Ati);	

	//auto A = douban::make_mat_csr(m, n, &Av, &Ap, &Ai);
	//auto At = douban::make_mat_csr(m_At, n_At, &Atv, &Atp, &Ati);

	//douban::linalg::svd_tr(A, U, S, V, p, 1.e-7, display);
	
	fclose(fp);
	MPI_Finalize();

	return 0;

}
