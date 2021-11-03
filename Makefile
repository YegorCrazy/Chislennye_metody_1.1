.PHONY: all clean again

all: main

main.o: main.c
	gcc -Wall -g -fsanitize=address -fsanitize=undefined -fsanitize=leak -c main.c -o main.o -lm

matrix.o: matrix.c
	gcc -Wall -g -fsanitize=address -fsanitize=undefined -fsanitize=leak -c matrix.c -o matrix.o -lm

main: main.o matrix.o
	gcc -Wall -g -fsanitize=address -fsanitize=undefined -fsanitize=leak matrix.o main.o -o main -lm

clean:
	rm -rf main *.o

again:
	make clean
	make all
