all: compila

compila: main.o 
	g++ main.o -o prog -Wl,--stack,2097152000

main.o: main.cpp
	g++ -c main.cpp

clean:
	rm main.o prog
