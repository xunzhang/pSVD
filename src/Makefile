#
# Makefile for svd test code
#

SUFFIX = .cpp
CC = mpic++
LIB = -L/home/wuhong_intern/local/lib
OPTS = -std=c++0x -std=gnu++0x -fopenmp -Wall -O3
FLAGS = -I/usr/include/mysql -I/home/wuhong_intern/local/include -llapack -lgfortran -lblas -lcblas -lmpi
SOURCE = main_svd_tr.o

all: main_svd_tr

main_svd_tr: $(SOURCE)
	rm -f main_svd_tr
	$(CC) -o main_svd_tr $(OPTS) $(FLAGS) $(LIBS) $(SOURCE)

$(SUFFIX).o:
	$(CC) $(OPTS) $(FLAGS) $(LIB) -c $*$(SUFFIX)

clean:
	/bin/rm -f main_svd_tr *.o
