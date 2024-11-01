Compile=g++ -g -std=c++20 -O3 # improve run time e.g. 10s to 4s
ideal_gas: main.o ideal_gas.o
	$(Compile) -o ideal_gas main.o ideal_gas.o
main.o: main.cpp ideal_gas.h
	$(Compile) -o main.o -c main.cpp
ideal_gas.o: ideal_gas.cpp ideal_gas.h
	$(Compile) -o ideal_gas.o -c ideal_gas.cpp
clean:
	rm -f main.o ideal_gas.o