/* test_vec_op.cpp */

#include "vec_op.hpp"
#include <vector>
#include <iostream>

int main(int argc, char *argv[]) 
{
	std::vector<double> a;
	a = vec_random<double>((size_t)2);
	for(std::vector<double>::size_type ix = 0; ix != a.size(); ++ix)
		std::cout << a[ix] << std::endl;
	std::cout << vec_unit(a) << std::endl;
	for(std::vector<double>::size_type ix = 0; ix != a.size(); ++ix)
		std::cout << a[ix] << std::endl;
	return 0;
}
