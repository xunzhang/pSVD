/* matrix_trans.cxx */

#include <map>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

struct value {
	int col;
	int row;
	double val;
};


typedef struct value Val;


int compare(const Val &a, const Val &b) {
	return a.col < b.col;
}


int main()
{
	int m, n, nnz;
	FILE *fp;
	/*	
	fp = fopen("pwtk.bin", "rb");
	fread(&m, sizeof(int), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	fread(&nnz, sizeof(int), 1, fp);
	
	std::vector<int> ar(nnz);
	std::vector<int> ac(nnz);
	std::vector<double> av(nnz);
	
	for(int i = 0; i < nnz; ++i) {
		fread(&(ar[i]), sizeof(int), 1, fp);
		fread(&(ac[i]), sizeof(int), 1, fp);
		fread(&(av[i]), sizeof(double), 1, fp);
	}
	*/
	
	m = 4;
	n = 3;
	nnz = 8;
	
	std::vector<int> ar;
	std::vector<int> ac;
	std::vector<double> av;
	
	
	ar.push_back(0);
	ar.push_back(0);
	ar.push_back(0);
	ar.push_back(1);
	ar.push_back(1);
	ar.push_back(2);
	ar.push_back(3);
	ar.push_back(3);
	ac.push_back(0);
	ac.push_back(1);
	ac.push_back(2);
	ac.push_back(0);
	ac.push_back(1);
	ac.push_back(2);
	ac.push_back(1);
	ac.push_back(2);
	
	av.push_back(1);
	av.push_back(2);
	av.push_back(3);
	av.push_back(4);
	av.push_back(5);
	av.push_back(6);
	av.push_back(7);
	av.push_back(8);
	
	std::vector<Val> obj;
	
	Val tmp;
	
	for(int i = 0; i < nnz; ++i) {
		tmp.col = ac[i];
		tmp.row = ar[i];
		tmp.val = av[i];
		obj.push_back(tmp);
	}
	for(int i = 0; i < nnz; ++i) 
		std::cout << obj[i].col << "\t" << obj[i].row << "\t" << obj[i].val << std::endl;;
	
	//qsort(obj.begin(), nnz, sizeof(Val), compare);
	std::sort(obj.begin(), obj.end(), compare);

	for(int i = 0; i < nnz; ++i) 
		std::cout << obj[i].col << "\t" << obj[i].row << "\t" << obj[i].val << std::endl;;
	
	//fclose(fp);
	return 0;
}
