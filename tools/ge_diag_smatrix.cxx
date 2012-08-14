/* ge_random_smatrix */
/* suppose row storage */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

double randdoub() {
	return 10 * (rand() / (double) RAND_MAX);
}

int main(int argc, char *argv[]) 
{
	int i, j;
	int m = 1000000, n = 1000000, nnz = 0;
	int row = 0, col;
	double val;
	FILE *fp = fopen("matrix_100wx100w.bin", "wb");

	nnz = 3 * m - 2;

	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	fwrite(&nnz, sizeof(int), 1, fp);

	for(i = 0; i < 2; ++i) {
		col = i;
		val = randdoub();
		fwrite(&row, sizeof(int), 1, fp);
		fwrite(&col, sizeof(int), 1, fp);
		fwrite(&val, sizeof(double), 1, fp);
	}
	for(i = 1; i < m - 1; ++i) {
		row = i;
		for(j = i - 1; j < i + 2; ++j) {
			col = j;
			val = randdoub();
			fwrite(&row, sizeof(int), 1, fp);
			fwrite(&col, sizeof(int), 1, fp);
			fwrite(&val, sizeof(double), 1, fp);
		}
	}
	row = m - 1;
	for(i = m - 2; i < m; ++i) {
		col = i;
		val = randdoub();
		fwrite(&row, sizeof(int), 1, fp);
		fwrite(&col, sizeof(int), 1, fp);
		fwrite(&val, sizeof(double), 1, fp);
	}

	fclose(fp);

	/* 
	 * print the random sparsity matrix.
	 */
	/*
	FILE *w_fp = fopen("matrix_1wx1w.bin", "rb");
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
