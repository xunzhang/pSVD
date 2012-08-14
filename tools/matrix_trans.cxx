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
	
	std::vector<Val> obj;
	
	Val tmp;
	
	for(int i = 0; i < nnz; ++i) {
		tmp.col = ac[i];
		tmp.row = ar[i];
		tmp.val = av[i];
		obj.push_back(tmp);
	}
	/*
	for(int i = 0; i < nnz; ++i) 
		std::cout << obj[i].col << "\t" << obj[i].row << "\t" << obj[i].val << std::endl;;
	*/
	//qsort(obj.begin(), nnz, sizeof(Val), compare);
	std::stable_sort(obj.begin(), obj.end(), compare);
	
	FILE *w_fp;
	w_fp = fopen("pwtk_trans.bin", "wb");
	fwrite(&m, sizeof(int), 1, w_fp);
	fwrite(&n, sizeof(int), 1, w_fp);
	fwrite(&nnz, sizeof(int), 1, w_fp);

	int tmp1, tmp2;
	double tmp3;
	for(int i = 0; i < nnz; ++i) {
		tmp1 = obj[i].col;
		tmp2 = obj[i].row;
		tmp3 = obj[i].val;
		std::cout << tmp1 << " " << tmp2 << std::endl;
		fwrite(&tmp1, sizeof(int), 1, w_fp);
		fwrite(&tmp2, sizeof(int), 1, w_fp);
		fwrite(&tmp3, sizeof(double), 1, w_fp);
	}

/*
	for(int i = 0; i < nnz; ++i) 
		std::cout << obj[i].col << "\t" << obj[i].row << "\t" << obj[i].val << std::endl;;
*/	
	fclose(fp);
	fclose(w_fp);
	return 0;
}
