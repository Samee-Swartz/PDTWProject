all: PDTW.exe

PDTW.exe: PDTW.o
	g++ -std=c++11 -Wall -g -o PDTW PDTW.o

PDTW.o: PDTW.cpp
	g++ -std=c++11 -Wall -c -g -o PDTW.o PDTW.cpp

clean:
	rm -f PDTW PDTW.o
