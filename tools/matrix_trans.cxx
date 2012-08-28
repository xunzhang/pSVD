/* matrix_trans.cxx */

#include <map>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <fstream>

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
	
	fp = fopen("db_data.txt", "r");
	fscanf(fp, "%d", &m);
	fscanf(fp, "%d", &n);
	fscanf(fp, "%d", &nnz);

	std::vector<int> ar(nnz);
	std::vector<int> ac(nnz);
	std::vector<double> av(nnz);
	
	for(int i = 0; i < nnz; ++i) {
		fscanf(fp, "%d", &ar[i]);
		fscanf(fp, "%d", &ac[i]);
		fscanf(fp, "%lf", &av[i]);
	}
	
	std::vector<Val> obj;
	
	Val tmp;
	
	for(int i = 0; i < nnz; ++i) {
		tmp.col = ac[i];
		tmp.row = ar[i];
		tmp.val = av[i];
		obj.push_back(tmp);
	}
	
	std::stable_sort(obj.begin(), obj.end(), compare);
	
	FILE *w_fp;
	w_fp = fopen("db_data_trans.txt", "w+");
	fprintf(w_fp, "%d ", m);
	fprintf(w_fp, "%d ", n);
	fprintf(w_fp, "%d\n", nnz);

	int tmp1, tmp2;
	double tmp3;
	for(int i = 0; i < nnz; ++i) {
		tmp1 = obj[i].col;
		tmp2 = obj[i].row;
		tmp3 = obj[i].val;
		fprintf(w_fp, "%d ", tmp1);
		fprintf(w_fp, "%d ", tmp2);
		fprintf(w_fp, "%lf\n", tmp3);
	}

	fclose(fp);
	fclose(w_fp);
	return 0;
}
