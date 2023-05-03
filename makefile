CXX = g++
VARNAME = value
CXXFLAGS = -Wall -g -Ofast -fopenmp

main: main.o NumMethods.o finitenuclei.o
	$(CXX) $(CXXFLAGS) -o main main.o NumMethods.o finitenuclei.o
	
NumMethods.o: NumMethods.hpp Conversions.hpp
	$(CXX) $(CXXFLAGS) -c NumMethods.cpp Conversions.cpp

Conversions.o: Conversions.hpp
	$(CXX) $(CXXFLAGS) -c Conversions.cpp

finitenuclei.o: finitenuclei.hpp NumMethods.hpp
	$(CXX) $(CXXFLAGS) -c finitenuclei.cpp NumMethods.cpp

clean:
	rm -f main main.o NumMethods.o finitenuclei.o Conversions.o