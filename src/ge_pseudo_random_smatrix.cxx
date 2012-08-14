/* ge_random_smatrix */
/* suppose row storage */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

int randint(int l, int u) {
	srand((unsigned)time(NULL));
	return l + rand() % (u - l + 1);
}

double randdoub() {
	return 10 * (rand() / (double) RAND_MAX);
}

void ShufKnuth(int nnz, int size, int *res_arry) {
	int i, j = 0;
	srand((unsigned int)time(NULL));
	for(i = 0; i < size; ++i) {
		if(rand() % (size - i) < nnz) {
			res_arry[j] = i;
			j++;
			nnz--;
		}
	}
}

int main(int argc, char *argv[]) 
{
	int i, j;
	int m = 10000, n = 10000, nnz = 0;
	double sparsity = 0.01;
	int *p, row, col, size;
	double val;
	FILE *fp = fopen("matrix_1wx1w.bin", "wb");
	
	nnz = m * n * sparsity;
	assert(nnz > m); // satisfy "sparsity * m > 1"
	
	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	fwrite(&nnz, sizeof(int), 1, fp);
	
	p = (int *)malloc(nnz * sizeof(int));
	size = nnz / m;
	for(i = 0; i < m - 1; ++i) {
		row = i;
		ShufKnuth(size, n, p);
		for(j = 0; j < size; ++j) {
			col = p[j];
			val = randdoub();
			fwrite(&row, sizeof(int), 1, fp);
			fwrite(&col, sizeof(int), 1, fp);
			fwrite(&val, sizeof(double), 1, fp);
		}
	}
	size = nnz - size * (m - 1);
	ShufKnuth(size, n, p);
	for(j = 0; j < size; ++j) {
		col = p[j];
		val = randdoub();
		fwrite(&row, sizeof(int), 1, fp);
		fwrite(&col, sizeof(int), 1, fp);
		fwrite(&val, sizeof(double), 1, fp);

	}	

	free(p);	
	fclose(fp);

	/* 
	 * print the random sparsity matrix.
	 */
/*
	FILE *w_fp = fopen("matrix_1kx1k.bin", "rb");
	int rint;
	double rdouble;

	fread(&rint, sizeof(int), 1, w_fp);
	printf("m is %d\t", rint);
	fread(&rint, sizeof(int), 1, w_fp);
	printf("n is %d\t", rint);
	fread(&rint, sizeof(int), 1, w_fp);
	printf("nnz is %d\n", rint);

	for(i = 0; i < nnz; ++i) {
		fread(&rint, sizeof(int), 1, w_fp);
		printf("row:%d\t", rint);
		fread(&rint, sizeof(int), 1, w_fp);
		printf("col:%d\t", rint);
		fread(&rdouble, sizeof(double), 1, w_fp);
		printf("val:%lf\n", rdouble);
	}
	fclose(w_fp);
*/

	return 0;
}
