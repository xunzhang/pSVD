#include <iostream>
#include <vector>
#include "mat_csr.hpp"
#include "mat_container.hpp"

int main(int argc, char* argv[]) 
{	
	size_t m = 6, n = 5, nnz = 10, k = 5, p = 4;
	
	int Ai_vec[10] = {0, 3, 1, 2, 4, 0, 1, 2, 3, 4}; // row_ptr size is nnz
	int Ap_vec[7] = {0, 2, 3, 5, 7, 9, 10}; // col_ind size is m + 1
	double Av_vec[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}; // val size is nnz
	
	std::vector<double> Av(Av_vec, Av_vec + nnz);
	std::vector<int> Ap(Ap_vec, Ap_vec + m + 1), Ai(Ai_vec, Ai_vec + nnz)   ;

	mat_csr A(m, n, nnz, Ai, Ap, Av);
	
	std::vector<double> U_tmp(m);
	std::vector< std::vector<double> > U(k, U_tmp); 
	/* mat_container<double> U(m, k) */
	mat_container<double> UU(m, k);
	
	std::vector<double> S(k);
	std::vector<double> S_tmp(n);
	std::vector< std::vector<double> > V(k, S_tmp); 
	
	std::cout << A.get_val_size() << std::endl;
	std::cout << U_tmp.size() << std::endl;
	std::cout << U.size() << std::endl;
	
	std::cout << UU.dim0() << "	" << UU.dim1() << std::endl;
	//svd_tr(A, U, S, V, p, 1.e-7);
	
	return 0;
}
