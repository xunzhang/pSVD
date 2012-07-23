/* mat_csr.hpp */
#ifndef MAT_CSR_HPP
#define MAT_CSR_HPP

#include <iostream>
#include <vector>
#include <cassert>

class mat_csr {
public:
	mat_csr(const size_t &m, const size_t &n, const size_t &num,
	        const std::vector<int> &r_p, 
		const std::vector<int> &c_i, 
		const std::vector<double> &v) : 
		dim0(m), 
		dim1(n), 
		nnz(num), 
		row_ptr(r_p), col_ind(c_i), val(v) {
			assert(row_ptr.size() == nnz);
			assert(col_ind.size() == (m + 1));
			assert(val.size() == nnz);
		}
/*
	mat_csr(const mat_csr &rhs) : 
		dim1(rhs.dim1), 
		dim2(rhs,dim2), 
		nnz(rhs.nnz), 
		row_ptr(&rhs.row_ptr), 
		col_ind(&rhs.col_ind),
		val(&rhs.val) {}
*/
	size_t get_row_ptr_size() { return row_ptr.size(); }
	size_t get_col_ind_size() { return col_ind.size(); }
	size_t get_val_size() { return val.size(); }


private:
	size_t dim0, dim1, nnz;
	std::vector<int> row_ptr; 
	std::vector<int> col_ind;
	std::vector<double> val;
};

#endif
