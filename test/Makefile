#
# Makefile for svd test code
#

SUFFIX = .cpp
CC = mpic++
LIB = -L/home/wuhong_intern/svd/lib
OPTS = -std=c++0x -std=gnu++0x -fopenmp -Wall -O3
CFLAG = -I/home/wuhong_intern/svd/include -I/usr/include/mysql -llapack -lgfortran -lblas -lcblas -lmpi
SOURCE = main_svd_tr.o

all: main_svd_tr

main_svd_tr: $(SOURCE)
	rm -f main_svd_tr
	$(CC) -o main_svd_tr $(OPTS) $(CFLAG) $(LIBS) $(SOURCE)

$(SUFFIX).o:
	$(CC) $(OPTS) $(CFLAG) $(LIB) -c $*$(SUFFIX)

clean:
	/bin/rm -f main_svd_tr *.o
