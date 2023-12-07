### possible array-based alternative to using splay tree for tour data structure?
**problem**: while the array-based representation is best for smaller problems, it is too slow for larger problems due to O(n) complexity per tour update.

**idea**: avoid long-distance segment reversals. How:
1. The initial tour heuristic should take care not to create situations that require many long-distance swaps to resolve. For this reason, the k-means method intuitively seems like a good fit.
2. On each pass of the LKH algo, skip swaps in a probabilistic manner: the probability of choosing to make the swap is a function of 2 things: the first is the number of swaps made so far in the pass, and the second is the sum of segment lengths to be reversed. Low distance swaps usually will not be skipped, and long distance swaps will be skipped with high probability unless a small number of swaps have been made previously in the current pass.

