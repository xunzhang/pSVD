#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{	
	int i = 10;
	MPI_Init(&argc, &argv);
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	printf("Hello world an d %d\n", i);
	MPI_Finalize();
}
