all: program count_kmers

findCommonKmers: findCommonKmers.cpp
	g++ -O3 $^ -o findCommonKmers

KMC/kmc_api/kmc_file.o:
	+make -C KMC -f makefile all

count_kmers: count_kmers_kmc.cpp KMC/kmc_api/kmc_file.o
	g++ -std=c++11 -Wall -I./KMC -I./KMC/kmc_counter ./KMC/kmc_api/*.o count_kmers_kmc.cpp -o count_kmers_kmc

program: program.cpp
	g++ -std=c++11 -DNDEBUG -O3 -I$(CPATH) $^ -o program -L$(LIBRARY_PATH) -lsdsl -ldivsufsort -ldivsufsort64

clean:
	rm -f findCommonKmers
	rm -f program
