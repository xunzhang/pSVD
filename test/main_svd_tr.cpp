/* run_svd_tr.cpp */

#include <map>
#include <mpi.h>
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <memory.h>
#include <douban/matvec.hpp>
#include <douban/option_parser.hpp>
#include <douban/linalg/gemv_ax.hpp>
#include <douban/matvec/mat_container.hpp>
#include <douban/matvec/vec_kahan_gemv.hpp>
#include <douban/mpi/linalg/mpi_svd_tr.hpp>

int main(int argc, char *argv[]) {

  int rank, nprocs;
  
  /* init MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  FILE *fp, *fp2;
  int m, n, m_At, n_At, nnz, k = 200, p = 200, display = 10;
  int len, begin, end, s_indx, t_s_indx;
  std::map<int, int> word_kind;
  
  std::string a_ascii_file = "../data/db_data.txt";
  std::string at_ascii_file = "../data/db_data_trans.txt";
	
  douban::option_parser cmd_parser;
  cmd_parser
      .add_help()
      .add_value_option("a-ascii-file", &a_ascii_file, "FILE::A in ascii format")
      .add_value_option("at-ascii-file", &at_ascii_file, "FILE::At in acsii format")
      ;
  cmd_parser.parse_all(argc, argv);
  
  //fp = fopen("../data/data.txt", "r");
  fp = fopen(a_ascii_file.c_str(), "r");
  fscanf(fp, "%d", &m);  fscanf(fp, "%d", &n);  fscanf(fp, "%d", &nnz);  fgetc(fp);
        
  //fp2 = fopen("../data/data_trans.txt", "r");
  fp2 = fopen(at_ascii_file.c_str(), "r");
  fscanf(fp2, "%d", &m_At);  fscanf(fp2, "%d", &n_At);  fscanf(fp2, "%d", &nnz);  fgetc(fp2);
	
  douban::mat_container<double> U(m, k);
  douban::mat_container<double> V(n, k);
  douban::vec_container<double> S(k);

  std::vector<size_t> ar, ac, Ap, Ai, atr, atc, Atp, Ati;
  std::vector<double> av, atv, Av, Atv;

  len = nnz / nprocs;	
  begin = rank * len;
  end = (rank + 1) * len;
  
  if(rank == nprocs - 1) {
    end = nnz;
    len = nnz - begin;
  }

  ar.resize(len);  ac.resize(len);  av.resize(len);
  atr.resize(len);  atc.resize(len);  atv.resize(len);
	
  char ch = '0';	
  for(int i = 0; i < len * rank; ++i) {
    ch = fgetc(fp);
    while(ch != '\n') { ch = fgetc(fp); }
  }	
  for(size_t i = 0; i < (size_t)len; ++i) {
    fscanf(fp, "%ld", &(ar[i]));
    fscanf(fp, "%ld", &(ac[i]));
    fscanf(fp, "%lf", &(av[i]));
  }
	
  
  char ch2 = '0';
  for(int i = 0; i < len * rank; ++i) {
    ch2 = fgetc(fp2);
    while(ch2 != '\n') { ch2 = fgetc(fp2); }
  }
  for(size_t i = 0; i < len; ++i) {
    fscanf(fp2, "%ld", &(atr[i]));
    fscanf(fp2, "%ld", &(atc[i]));
    fscanf(fp2, "%lf", &(atv[i]));
  }

  m = 0;
  m_At = 0;
  for(size_t i = 0; i < (size_t)len; ++i) { 
    if(ar[i] != ar[i+1])
      m++;
    if(atr[i] != atr[i+1])
      m_At++;
  }
/*
  m = ar[len - 1] - ar[0] + 1;
  m_At = atr[len - 1] - atr[0] + 1;
*/
  s_indx = ar[0];
  t_s_indx = atr[0];
  if(rank != 0) {
    for(size_t i = 0; i < (size_t)len; ++i) ar[i] -= s_indx;
    for(size_t i = 0; i < (size_t)len; ++i) atr[i] -= t_s_indx;	
  }

  Ap.resize(m + 1);  Av.resize(len);  Ai.resize(len);
  Atp.resize(m_At + 1);  Atv.resize(len);  Ati.resize(len);
        
  // transfer coo to csr
  douban::coo_to_csr((size_t)m, (size_t)n, (size_t)len, av, ar, ac, Av, Ap, Ai);	
  douban::coo_to_csr((size_t)m_At, (size_t)n_At, (size_t)len, atv, atr, atc, Atv, Atp, Ati);	
  
  // init A and A'
  auto A = douban::make_mat_csr((size_t)m, (size_t)n, &Av, &Ap, &Ai);
  auto At = douban::make_mat_csr((size_t)m_At, (size_t)n_At, &Atv, &Atp, &Ati);
  
  // svd fact
  douban::linalg::svd_tr(A, At, U, S, V, p, 1.e-7, display, s_indx, t_s_indx);
  
  //douban::linalg::svd_tr(douban::linalg::make_gemv_ax(&A), douban::linalg::make_gemv_ax(&At), U, S, V, p, 1.e-7, display, s_indx, t_s_indx);
  /*
  if(rank == 0) {
    std::cout << "svd_tr finished!" << std::endl;
    double cost = douban::linalg::svd_cost(A, U, S, V);
    std::cout << "cost of SVD is " << cost << std::endl; 
  }
  */
  
  // print Sigular-Values
  for(size_t i = 0; i < (size_t)k; i++)
    std::cout << "S(i) is" << S.get(i) << std::endl;
        	
  fclose(fp);
  fclose(fp2);
  MPI_Finalize();

  return 0;
}
