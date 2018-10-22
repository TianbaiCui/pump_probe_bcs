MKLROOT = /opt/local
CXXFLAGS = -std=gnu++14 -I$(MKLROOT)/include
CXX = g++


main: main.o
	$(CXX) $(CXXFLAGS) main.o -o main

main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) main.cpp

main_dense: main_dense.o
	$(CXX) $(CXXFLAGS) main_dense.o -o main_dense

main_dense.o : main_dense.cpp
	$(CXX) -c $(CXXFLAGS) main_dense.cpp

main_check: main_check.o
	$(CXX) $(CXXFLAGS) main_check.o -o main_check

main_check.o : main_check.cpp
	$(CXX) -c $(CXXFLAGS) main_check.cpp

main_A: main_A.o
	$(CXX) $(CXXFLAGS) main_A.o -o main_A

main_A.o : main_A.cpp
	$(CXX) -c $(CXXFLAGS) main_A.cpp

clean:
	rm -f *.o *~main *~main_A *~main_dense
