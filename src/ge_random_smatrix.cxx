/* ge_random_smatrix */
/* suppose row storage */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
  	int i;
	int m = 10, n = 10, nnz = 0;
	double sparsity = 0.1;
	int *p, row, col;
	double val;
	FILE *fp = fopen("matrix_demo1.bin", "wb");
	
	nnz = m * n * sparsity;
	
	fwrite(&m, sizeof(int), 1, fp);
	fwrite(&n, sizeof(int), 1, fp);
	fwrite(&nnz, sizeof(int), 1, fp);
	
  	p = (int *)malloc(nnz * sizeof(int));
	
	ShufKnuth(nnz, m * n, p);
	
	for(i = 0; i < nnz; ++i) {
		row = p[i] / n;
		col = p[i] % n;
		val = randdoub();
		fwrite(&row, sizeof(int), 1, fp);
		fwrite(&col, sizeof(int), 1, fp);
		fwrite(&val, sizeof(double), 1, fp);
	}
	
	free(p);	
	fclose(fp);
	printf("%lf\n", randdoub());
	/* 
	 * print the random sparsity matrix.
	 */
	
	/*	
	FILE *w_fp = fopen("matrix_demo1.bin", "rb");
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
