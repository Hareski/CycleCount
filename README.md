# Algorithm for Counting Cycles in C++

This is a general purpose algorithm for counting simple cycles of up to a fixed length included on an adjacency matrix passing through all or a specified vertex.
Working only on an undirected graph but will be updated soon.

Based on this article: Giscard, P., Kriege, N. & Wilson, R.C. A General Purpose Algorithm for Counting Simple Cycles and Simple Paths of Any Length. Algorithmica 81, 2716–2737 (2019). https://doi.org/10.1007/s00453-019-00552-1

Some others implementations of this algorithm:
 * Matlab: Pierre-Louis Giscard - https://www.mathworks.com/matlabcentral/fileexchange/60814-cyclecount-a-l0
 * Python: Morteza Milani - https://github.com/milani/cycleindex

## Some examples
### Counting simple cycle in the whole graph
```c++
int main() {
    std::vector<std::vector<int>> G = {{0, 1, 0, 0, 1},
                                       {1, 0, 1, 1, 1},
                                       {0, 1, 0, 1, 0},
                                       {0, 1, 1, 1, 1},
                                       {1, 1, 0, 1, 0}};
    const std::vector<double> result = cycleCount(G, 5, false);
    for (const auto c : result) {
        std::cout << std::setw(15) << std::left << round(c);
    }
}
```
Display the array `1  7  6  4  2` because the graph described by the adjacency matrix G contains 1 self-loop, 7 cycles of length 2 (just the edges) and respectively 6, 4 and 2 cycles of length 3, 4 and 5.

### Counting simple cycle passing trough a fixed vertex
It is possible to get count of all simple cycles passing trough a fixed vertex having as index 3:
```c++
    const std::vector<double> countFixedVertex = cycleCountFixedVertex(adjacencyMatrix, 3, adjacencyMatrix.size());
```

## Perspectives
It is possible to extend the use of eigenvalues to compute the trace of directed graph or use a better algorithm to do the matrix multiplication for a fixed vertex.

**THE ALGORITHM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED**