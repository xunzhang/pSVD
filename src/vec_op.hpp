/* vec_op.hpp */
#ifndef VEC_OP_HPP
#define VEC_OP_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>

template<typename T> T norm2(const std::vector<T> &vec) {
	T n = 0.0;  // is this right?
	for(typename std::vector<T>::size_type ix = 0; ix != vec.size(); ++ix) {
		n += vec[ix] * vec[ix];
	}	
	return (T) sqrt(n);
} 

template<typename T> void scale(std::vector<T> &vec, const T &coff) {
	for(typename std::vector<T>::size_type ix = 0; ix != vec.size(); ++ix) {
		vec[ix] *= coff;
	}
	return;
}

template<typename T> void vec_random(size_t n, std::vector<T> &res) {
	T rand_tmp;
	//std::vector<T> res(n, 0);
	srand((unsigned)time(NULL));
	for(size_t i = 0; i < n; ++i) {
		rand_tmp = (T)((T)rand()/(T)RAND_MAX);
		res.push_back(rand_tmp);
	}
	return;
}

template<typename T> T vec_unit(std::vector<T> &vec) {
	T norm = norm2(vec);
	scale(vec, (T)(1.0 / norm));
	return norm;
}

#endif
