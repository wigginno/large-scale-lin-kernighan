# large-scale-lin-kernighan

## Project description

TSP heuristic solver designed to be used on very large-scale instances like GalaxyTSP[^1] using geometric methods and local search.

## Progress

**Completed**
- Prototype initial tour generation using k-means clustering, runs in  O(nlogn), linear space
- Prototype Delaunay triangulation (for candidate set generation) using randomized incremental insertion, runs in O(nlogn), linear space
- Tested initial tour generation and Delaunay triangulation on 100 million city problem

**In-progress**
- Partition scheme for n > 100 million cities
- Transform partitions of geospatial positions to Cartesian coordinates
- Splay tree data structure to represent tour[^2]

[^1]: Drori, I., Kates, B., Sickinger, W., Kharkar, A., Dietrich, B., Shporer, A., & Udell, M. (2020). GalaxyTSP: A New Billion-Node Benchmark for TSP. 1st Workshop on Learning Meets Combinatorial Algorithms @ NeurIPS 2020, Vancouver, Canada. Retrieved from [https://www.cs.columbia.edu/~idrori/galaxytsp.pdf](https://www.cs.columbia.edu/~idrori/galaxytsp.pdf)

[^2]: Fredman, M., Johnson, D., Mcgeoch, L., & Ostheimer, G. (1995). Data structures for traveling salesmen. Journal of Algorithms, 18(3), 432-479. [https://doi.org/10.1006/jagm.1995.1018](https://doi.org/10.1006/jagm.1995.1018)
