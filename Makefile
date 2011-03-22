CC=g++

all: aligator_exe

aligator_exe: aligator.o common.o pathway_db.o
	$(CC) -O3 -g -Wall  -o aligator_exe aligator.o common.o pathway_db.o libboost_program_options-gcc43-mt-1_35.a

common.o: common.cpp common.hpp
	$(CC) -O3 -g -Wall -c -I. common.cpp

aligator.o: aligator.cpp aligator.hpp
	$(CC) -O3 -g -Wall -c -I. aligator.cpp

pathway_db.o: pathway_db.cpp pathway_db.hpp
	$(CC) -O3 -g -Wall -Wno-write-strings -c -I. pathway_db.cpp 

clean:
	rm *.o aligator_exe
