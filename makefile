CXX = g++-10
CC = gcc-10
VARNAME = value
CXXFLAGS = -Wall -g -O3 -fopenmp -std=c++17 -fPIC

main: main.o NumMethods.o finitenuclei.o MCMC.o infinitematter.o Conversions.o minpack.o
	$(CXX) $(CXXFLAGS) -o main main.o NumMethods.o finitenuclei.o MCMC.o infinitematter.o Conversions.o minpack.o

# Target for the shared library
cpylibrary.so: main.o NumMethods.o finitenuclei.o MCMC.o infinitematter.o Conversions.o minpack.o
	$(CXX) $(CXXFLAGS) -dynamiclib -o cpylibrary.so main.o NumMethods.o finitenuclei.o MCMC.o infinitematter.o Conversions.o minpack.o

NumMethods.o: NumMethods.hpp Conversions.hpp
	$(CXX) $(CXXFLAGS) -c NumMethods.cpp Conversions.cpp

Conversions.o: Conversions.hpp
	$(CXX) $(CXXFLAGS) -c Conversions.cpp

finitenuclei.o: finitenuclei.hpp NumMethods.hpp
	$(CXX) $(CXXFLAGS) -c finitenuclei.cpp NumMethods.cpp

MCMC.o: MCMC.hpp NumMethods.hpp finitenuclei.hpp infinitematter.hpp Conversions.hpp
	$(CXX) $(CXXFLAGS) -c MCMC.cpp NumMethods.cpp finitenuclei.cpp infinitematter.cpp Conversions.cpp

infinitematter.o: infinitematter.hpp NumMethods.hpp
	$(CXX) $(CXXFLAGS) -c infinitematter.cpp NumMethods.cpp

minpack.o: minpack.hpp
	$(CXX) $(CXXFLAGS) -c minpack.cpp
	
clean:
	rm -f main main.o NumMethods.o finitenuclei.o Conversions.o MCMC.o infinitematter.o minpack.o cpylibrary.so