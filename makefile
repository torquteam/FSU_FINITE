CXX = g++
VARNAME = value
CXXFLAGS = -Wall -g -Ofast -fopenmp

main: main.o NumMethods.o finitenuclei.o MCMC.o infinitematter.o Conversions.o
	$(CXX) $(CXXFLAGS) -o main main.o NumMethods.o finitenuclei.o MCMC.o infinitematter.o Conversions.o
	
NumMethods.o: NumMethods.hpp Conversions.hpp
	$(CXX) $(CXXFLAGS) -c NumMethods.cpp Conversions.cpp

Conversions.o: Conversions.hpp
	$(CXX) $(CXXFLAGS) -c Conversions.cpp

finitenuclei.o: finitenuclei.hpp NumMethods.hpp
	$(CXX) $(CXXFLAGS) -c finitenuclei.cpp NumMethods.cpp

MCMC.o: MCMC.hpp NumMethods.hpp finitenuclei.hpp
	$(CXX) $(CXXFLAGS) -c MCMC.cpp NumMethods.cpp finitenuclei.cpp

infinitematter.0: infinitematter.hpp NumMethods.hpp
	$(CXX) $(CXXFLAGS) -c infinitematter.cpp NumMethods.cpp

clean:
	rm -f main main.o NumMethods.o finitenuclei.o Conversions.o MCMC.o infinitematter.o