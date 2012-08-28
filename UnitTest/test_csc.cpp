/* run_svd_tr.cpp */

/* Usage: mpic++ -I/home/wuhong_intern/local/include -L/home/wuhong_intern/local/lib -std=c++0x -std=gnu++0x -lmpi -fopenmp -I/usr/include/mysql      -Wall -O0 main.cpp -llapack -lgfortran -lblas -lcblas
*/

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

int main() 
{
	std::vector<int> ar, ac, Ap(3), Ai(3);
	std::vector<double> av, Av(3);
	
	
	ar.push_back(0);
	ar.push_back(1);
	ar.push_back(0);
	
	ac.push_back(0);
	ac.push_back(0);
	ac.push_back(1);
	
	av.push_back(1);
	av.push_back(2);
	av.push_back(3);
	
	douban::coo_to_csr(2, 4, 3, av, ar, ac, Av, Ap, Ai);	
	
	for(int i = 0; i < 3; ++i)
		std::cout << "Ap " << Ap[i] << std::endl;
	for(int i = 0; i < 3; ++i) {
		std::cout << "Ai " << Ai[i] << std::endl;
		std::cout << "Av" << Av[i] << std::endl;
	}
	return 0;
}
