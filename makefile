all: PDTW

PDTW: PDTW.o
	g++ -std=c++11 -Wall -o PDTW PDTW.o

PDTW.o: PDTW.cpp
	g++ -std=c++11 -Wall -c -o PDTW.o PDTW.cpp

clean:
	rm -f PDTW PDTW.o
