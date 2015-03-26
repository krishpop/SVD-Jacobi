CC = gcc

all: jac clean

jac: jacobi.o
	${CC} -o $@ $^

clean:
	rm -f *.o
