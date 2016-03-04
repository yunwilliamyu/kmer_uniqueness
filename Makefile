IDIR=.
CC=g++
LIB=/home/ywy/bin/lib
#CFLAGS=-I${IDIR} -I/usr/local/includes -std=c++11 -L${LIB} -Wall -O3 -DNDEBUG
CFLAGS=-I${IDIR} -I/home/ywy/bin/include -std=c++11 -pipe -L${LIB} -Wall -O3 -fopenmp \
#CFLAGS=-I${IDIR} -I/home/ywy/bin/include -std=c++11 -pipe -L${LIB} -Wall -O3 -DNDEBUG -fopenmp \
#CFLAGS=-I${IDIR} -I/home/ywy/bin/include -std=c++11 -pipe -march=native -L${LIB} -Wall -O3 -DNDEBUG -fopenmp \
# -pg \
	   #-ftree-vectorizer-verbose=2

OBJS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: ${OBJS} ref_distance_variable
	echo "All made."

ref_distance_variable: ref_distance_variable.o
	${CC} -o $@ ${CFLAGS} $@.o

%.o: %.cpp 
	${CC} ${CFLAGS} -c -o $@ $<

clean:
	rm -f ${OBJS} ref_distance_variable
	@echo "All cleaned up!"
