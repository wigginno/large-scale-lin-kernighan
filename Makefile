CXX = g++
CXXFLAGS = -std=c++17 -Wall
DEBUGFLAGS = -g
RELEASEFLAGS = -O3
TARGETS = kmeans_tsp gen_tsp delaunay

all: $(TARGETS)

debug: CXXFLAGS += $(DEBUGFLAGS)
debug: $(TARGETS)

release: CXXFLAGS += $(RELEASEFLAGS)
release: $(TARGETS)

kmeans_tsp: kmeans_tsp.o
	$(CXX) $(CXXFLAGS) -o kmeans_tsp kmeans_tsp.o

kmeans_tsp.o: kmeans_tsp.cpp kmeans_tsp.h
	$(CXX) $(CXXFLAGS) -c kmeans_tsp.cpp

delaunay: delaunay.o
	$(CXX) $(CXXFLAGS) -o delaunay delaunay.o

delaunay.o: delaunay.cpp delaunay.h
	$(CXX) $(CXXFLAGS) -c delaunay.cpp

gen_tsp: gen_tsp.o
	$(CXX) $(CXXFLAGS) -o gen_tsp gen_tsp.o

gen_tsp.o: gen_tsp.cpp
	$(CXX) $(CXXFLAGS) -c gen_tsp.cpp

clean:
	rm -f *.o $(TARGETS)
