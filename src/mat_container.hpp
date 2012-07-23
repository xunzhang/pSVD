/* mat_container.hpp */
#ifndef MAT_CONTAINER_HPP
#define MAT_CONTAINER_HPP

#include <vector>

template <class Type> class mat_container {
public:
	mat_container(size_t d0, size_t d1) {
		vec_dim0 = d0;
		vec_dim1 = d1;
		inter_vec.resize(vec_dim0);
		self.resize(vec_dim1, inter_vec);
	}
	size_t dim0() { return vec_dim0; }
	size_t dim1() { return vec_dim1; }
private:
	std::vector<Type> inter_vec;
	std::vector< std::vector<Type> > self;
	size_t vec_dim0, vec_dim1;
};

#endif
