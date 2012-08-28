#include <stdio.h>
#include <iostream>
int main(void)
{
	FILE *fp = fopen("pwtk_trans.bin", "rb");
	int i;
	double d;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&d, sizeof(double), 1, fp);
	std::cout << d << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&d, sizeof(double), 1, fp);
	std::cout << d << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&i, sizeof(int), 1, fp);
	std::cout << i << std::endl;
	fread(&d, sizeof(double), 1, fp);
	std::cout << d << std::endl;
}
