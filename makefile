main.x: main.c
	gcc -o main.x main.c -lm

build: main.x

run: 
	./main.x