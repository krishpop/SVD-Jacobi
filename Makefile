CC = gcc

all: hw2 clean

hw2: hw2.o
	${CC} -o $@ $^

clean:
	rm -f *.o
