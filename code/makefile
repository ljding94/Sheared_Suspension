Compile=g++ -g -std=c++20 -O3 # improve run time e.g. 10s to 4s
suspension: main.o suspension.o
	$(Compile) -o suspension main.o suspension.o
main.o: main.cpp suspension.h
	$(Compile) -o main.o -c main.cpp
suspension.o: suspension.cpp suspension.h
	$(Compile) -o suspension.o -c suspension.cpp
clean:
	rm -f main.o suspension.o