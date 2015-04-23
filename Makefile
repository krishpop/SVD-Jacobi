CC = gcc

all: jac clean

jac: jacobi.o
	${CC} -o $@ $^

jacobi.o: MAT_FUNC.h

clean:
	rm -f *.o
