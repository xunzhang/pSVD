/* test_mpi_func_call.cxx */

#include <stdio.h>
#include <mpi.h>

void mpi_func() {
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	printf("Hello world! This is process %d out of %d\n", rank, nprocs);
	if(rank == 0)
		printf("Some processes are more equal than others.");
	MPI_Finalize();
}

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	mpi_func();
	return 0;
}
