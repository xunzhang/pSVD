/* mat_container.hpp */
#ifndef MAT_CONTAINER_HPP
#define MAT_CONTAINER_HPP

#include <vector>

template <typename Type, bool row_major=true> class mat_container {
public:
	mat_container(size_t d0, size_t d1) {
		vec_dim0 = d0;
		vec_dim1 = d1;
		self.resize(vec_dim0 * vec_dim1);
	}

	size_t dim0() { return vec_dim0; }
	
	size_t dim1() { return vec_dim1; }
	
	inline void operator=(const Type val) {
		self.assign(vec_dim0 * vec_dim1, val);
	}
	
	std::vector<Type> col(size_t indx) {
		std::vector<Type> col;
		col.resize(vec_dim1);
		if(row_major) {
			for(std::vector<Type>::size_type ix = indx; ix != self.size(); ix += vec_dim1) { col.push_back(self[ix]); }
		} else {
			for(int i = 0, std::vector<Type>::size_type ix = indx * vec_dim0; i < vec_dim0; ++i, ++ix) { col.push_back(self[ix]); }
		}	
		return col;
	}

private:
	std::vector<Type> self;
	size_t vec_dim0, vec_dim1;
};

#endif
