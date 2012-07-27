/* mat_coo.hpp */
#ifndef MAT_COO_HPP
#define MAT_COO_HPP

#include <vector>

class mat_coo {
public:
	mat_coo(const size_t &m, const size_t &n, const size_t &num,
	 	const std::vector<size_t> &r_indx,
		const std::vector<size_t> &c_indx,
		const std::vector<double> &v) :
		dim0(m),
		dim1(n),
		nnz(num),
		row_indx(r_indx),
		col_indx(c_indx),
		val(v) {
			assert(row_indx.size() == nnz);
			assert(col_indx.size() == nnz);
			assert(val.size() == nnz);
		}
	mat_coo trans() {
		mat_coo At(dim1, dim0, nnz, col_indx, row_indx, val);
		return At;
	}

private:
	size_t dim0, dim1, nnz;
	std::vector<size_t> row_indx;
	std::vector<size_t> col_indx;
	std::vector<double> val;
};

#endif
