/* vec_op.hpp */
#ifndef VEC_OP_HPP
#define VEC_OP_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <stdlib.h>

template<typename T> double norm2(const std::vector<T> &vec) {
	double n = 0.0;
	for(std::vector<double>::size_type ix = 0; ix != vec.size(); ++ix) {
		n += vec[ix] * vec[ix];
	}	
	return sqrt(n);
} 

template<typename T> void scale(std::vector<T> &vec, double coff) {
	for(std::vector<double>::size_type ix = 0; ix != vec.size(); ++ix) {
		vec[ix] *= coff;
	}
	return;
}

template<typename T> std::vector<T> vec_random(size_t n) {
	T rand_tmp;
	std::vector<T> res;
	
	srand((unsigned)time(NULL));
	
	for(size_t i = 0; i < n; ++i) {
		rand_tmp = (T)((T)rand()/(T)RAND_MAX);
		res.push_back(rand_tmp);
	}
	return res;
}

template<typename T> double vec_unit(std::vector<T> &vec) {
	double norm = norm2<double>(vec);
	scale(vec, 1.0 / norm);
	return norm;
}

#endif
