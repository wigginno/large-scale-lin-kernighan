# large-scale-lin-kernighan

## Project description

A TSP heuristic solver designed for large-scale instances such as GalaxyTSP[^1].

## Progress

**Completed**
- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/wigginno/large-scale-lin-kernighan/blob/main/toy_lkh.ipynb)
 Interactive proof of concept for the LKH algorithm: improves tours to around 2% of optimal for small problems (up to 3,000 cities) **in milliseconds**. Note that an array based data structure will not be used in the real implementation, as it does not scale well for large problems.

**In-progress**
- Implement splay tree data structure to represent tour[^2], achieve logarithmic amortized complexity for tour update and lookup.
- Implement candidate set generation
- Implement initial tour generation

[^1]: Drori, I., Kates, B., Sickinger, W., Kharkar, A., Dietrich, B., Shporer, A., & Udell, M. (2020). GalaxyTSP: A New Billion-Node Benchmark for TSP. 1st Workshop on Learning Meets Combinatorial Algorithms @ NeurIPS 2020, Vancouver, Canada. Retrieved from [https://www.cs.columbia.edu/~idrori/galaxytsp.pdf](https://www.cs.columbia.edu/~idrori/galaxytsp.pdf)

[^2]: Fredman, M., Johnson, D., Mcgeoch, L., & Ostheimer, G. (1995). Data structures for traveling salesmen. Journal of Algorithms, 18(3), 432-479. [https://doi.org/10.1006/jagm.1995.1018](https://doi.org/10.1006/jagm.1995.1018)
