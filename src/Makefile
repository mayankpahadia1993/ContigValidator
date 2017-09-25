all: findCommonKmers program

findCommonKmers: findCommonKmers.cpp
	g++ -O3 $^ -o findCommonKmers

program: program.cpp
	g++ -std=c++11 -DNDEBUG -O3 -I$(CPATH) $^ -o program -L$(LIBRARY_PATH) -lsdsl -ldivsufsort -ldivsufsort64

clean:
	rm -f findCommonKmers
	rm -f program
