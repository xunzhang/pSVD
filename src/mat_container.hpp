/* mat_container.hpp */
#ifndef MAT_CONTAINER_HPP
#define MAT_CONTAINER_HPP

#include <vector>

template <typename T, bool row_major=true> class mat_container {
public:
	mat_container(size_t d0, size_t d1) {
		vec_dim0 = d0;
		vec_dim1 = d1;
		self.resize(vec_dim0 * vec_dim1);
	}
	
	mat_container(double *addr, size_t d0, size_t d1) {
		//self[0] = addr;
		vec_dim0 = d0;
		vec_dim1 = d1;
		self.resize(vec_dim0 * vec_dim1);
		self.assign(addr, addr + vec_dim0 * vec_dim1);
	}

	size_t dim0() { return vec_dim0; }
	
	size_t dim1() { return vec_dim1; }
	
	inline void operator=(const T &val) {
		self.assign(vec_dim0 * vec_dim1, val);
	}
	
	inline T operator()(size_t indx0, size_t indx1) {
		if(row_major)
			return self[indx0 * vec_dim1 + indx1];
		else
			return self[indx1 * vec_dim0 + indx0];
	}

	std::vector<T> col(size_t indx) {
		std::vector<T> col;
		col.resize(vec_dim1);
		if(row_major) {
			for(size_t ix = indx; ix != self.size(); ix += vec_dim1) { col.push_back(self[ix]); }
		} else {
			size_t ix = indx * vec_dim0;
			for(size_t i = 0; i < vec_dim0; i++, ix++) { col.push_back(self[ix]); }
		}	
		return col;
	}

private:
	std::vector<T> self;
	size_t vec_dim0, vec_dim1;
};

#endif
