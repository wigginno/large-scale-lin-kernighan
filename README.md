# large-scale-lin-kernighan

## Project description

A TSP heuristic solver designed for large-scale instances such as GalaxyTSP[^1], using geometric methods and local search.

## Progress

**Completed**
- Prototype initial tour generation using k-means clustering
- Prototype Delaunay triangulation (for candidate set generation) using randomized incremental insertion
- Tested initial tour generation and Delaunay triangulation on 100 million city problem

**In-progress**
- Transform partitions of geospatial positions to Cartesian coordinates
- Implement splay tree data structure to represent tour[^2]
- Proper source code organization - header files for declarations and cpp files for implementation
- Reimplement Delaunay triangulation w/ streaming technique to generate 1bn+ triangle meshes with low memory footprint

## Setup

```zsh
git clone https://github.com/wigginno/large-scale-lin-kernighan.git
cd large-scale-lin-kernighan
make release
```

## Usage

### Problem file generation
```zsh
./gen_tsp [n] [output_file]
```

### Initial tour generation
```zsh
./kmeans_tsp [problem file]
```

### Delaunay triangulation
```zsh
./delaunay [problem file]
```

## Algorithmic efficiency notes (in progress)
#### Initial tour generation
Runs in O(nlogn) time and O(n) space.
#### Generating next city candidates for tour improvement (Delaunay triangulation and post-processing)
Runs in O(nlogn) time and O(n) space.
#### Tour improvement (local search heuristic, preset number of iterations)
O(nlogn) expected time and O(n) space.

## Experimental results
in progress...

[^1]: Drori, I., Kates, B., Sickinger, W., Kharkar, A., Dietrich, B., Shporer, A., & Udell, M. (2020). GalaxyTSP: A New Billion-Node Benchmark for TSP. 1st Workshop on Learning Meets Combinatorial Algorithms @ NeurIPS 2020, Vancouver, Canada. Retrieved from [https://www.cs.columbia.edu/~idrori/galaxytsp.pdf](https://www.cs.columbia.edu/~idrori/galaxytsp.pdf)

[^2]: Fredman, M., Johnson, D., Mcgeoch, L., & Ostheimer, G. (1995). Data structures for traveling salesmen. Journal of Algorithms, 18(3), 432-479. [https://doi.org/10.1006/jagm.1995.1018](https://doi.org/10.1006/jagm.1995.1018)
