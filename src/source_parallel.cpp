#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <memory.h>

int main(int argc, char *argv[]){
	int numprocs, myid;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	int M, N, nnz;
	int *I, *J, *row_start;
	double *val;
	FILE *fp, *fp_in, *fp_out, *fp_test;
	double t_start = 0.0, t_end = 0.0;
	int i, j, begin, end, len;
	double *x, *b, *test, *result; 
	if(myid == 0){t_start = MPI_Wtime();}
	
	fp = fopen("pwtk.bin", "rb");
	
	fread((char*)&M, sizeof(int),1,fp);
	fread((char*)&N, sizeof(int),1,fp);
	fread((char*)&nnz, sizeof(int),1,fp);
	
	printf("\n");
	printf("Matrix file general info:\n");
	printf("Dimension size: %d\tx %d\n",M,N);
	printf("Number of nonzeros: %d\n",nnz);
	
	len = nnz/numprocs;
	begin = myid * len;
	end = (myid+1) * len;
	if(myid == numprocs - 1){
		end = nnz;
		len = nnz - begin;
	}
 
	I = malloc(len*sizeof(int));
	J = malloc(len*sizeof(int));
	val = malloc(len*sizeof(double));
	x = malloc(M*sizeof(double));
	b = malloc(M*sizeof(double));
	test = malloc(M*sizeof(double));
	memset(b,0,M*sizeof(double));
	if (myid == 0) {
		result = malloc(M*numprocs*sizeof(double));
	}	
	fseek(fp,(long)(begin*(sizeof(int)*2+sizeof(double))),SEEK_CUR);
	for(i = 0; i < len; i++){
		fread(&(I[i]), sizeof(int),1,fp);
		fread(&(J[i]), sizeof(int),1,fp);
		fread(&(val[i]), sizeof(double),1,fp);
	}
	
	fp_in = fopen("input.bin","rb");
	for(i = 0; i < M; i++){
		fread(&x[i],sizeof(double),1,fp_in);
	}
  
	for (i = 0; i < len; i++) {
		b[I[i]] += val[i] * x[J[i]];
	}
	
	MPI_Gather(b,M,MPI_DOUBLE,result,M,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
	if (myid == 0) {
		for (i =0; i < M; i++) {
			for (j = 1; j < numprocs; j++) {
				result[i] += result[j*M+i];
			}
		}
		if((fp_out = fopen("result.bin","w+")) == NULL){
			printf("Cannot open!\n");
			exit(0);
		}
		
		for(i = 0; i < M; i++){
			fwrite(&result[i],sizeof(double),1,fp_out);
		}
		t_end = MPI_Wtime();
		printf("Time is %.6f ms. \n",(t_end-t_start)*1000);
	}
	free(I);
	free(J);
	free(val);
	free(x);
	free(b);
	MPI_Finalize();
	return 0;
}

