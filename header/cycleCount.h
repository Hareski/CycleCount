#ifndef CYCLECOUNT_H
#define CYCLECOUNT_H

std::vector<std::vector<std::vector<double>>> pathsCount(std::vector<std::vector<double>> &adjacencyMatrix,
                                                             unsigned long length, bool directed);

std::vector<std::vector<std::vector<double>>> recursiveSubgraphsPaths(std::vector<std::vector<double>> adjacencyMatrix,
                                                                      unsigned long length,
                                                                      std::vector<int> subgraph,
                                                                      std::vector<bool> allowedVertex,
                                                                      std::vector<std::vector<std::vector<double>>> primes,
                                                                      std::vector<bool> neighbourhood,
                                                                      bool directed);

std::vector<std::vector<std::vector<double>>> primeCountPath(const std::vector<std::vector<double>> &adjacencyMatrix,
                                                             unsigned long length,
                                                             std::vector<int> subgraph, unsigned long neighboursNumber,
                                                             std::vector<std::vector<std::vector<double>>> primes);

std::vector<double> cycleCount(std::vector<std::vector<double>> &adjacencyMatrix,
                               unsigned long length, bool directed);

std::vector<double> recursiveSubgraphs(std::vector<std::vector<double>> adjacencyMatrix,
                                       unsigned long length,
                                       std::vector<int> subgraph,
                                       std::vector<bool> allowedVertex, std::vector<double> primes,
                                       std::vector<bool> neighbourhood,
                                       bool directed);

std::vector<double> primeCount(const std::vector<std::vector<double>> &adjacencyMatrix,
                               unsigned long length,
                               std::vector<int> subgraph,
                               unsigned long neighboursNumber,
                               std::vector<double> primes,
                               bool directed);

std::vector<double>
primeCountUndirected(const std::vector<std::vector<double>> &adjacencyMatrix, unsigned long length,
                     std::vector<int> subgraph, unsigned long neighboursNumber,
                     std::vector<double> primes);


std::vector<double>
primeCountDirected(const std::vector<std::vector<double>> &adjacencyMatrix, unsigned long length,
                   std::vector<int> subgraph, unsigned long neighboursNumber,
                   std::vector<double> primes);

std::vector<double> cycleCountFixedVertex(std::vector<std::vector<double>> &adjacencyMatrix,
                                          int vertex,
                                          unsigned long length);

std::vector<double> recursiveSubgraphsFixedVertex(std::vector<std::vector<double>> adjacencyMatrix,
                                                  int vertex,
                                                  unsigned long length,
                                                  std::vector<int> subgraph,
                                                  std::vector<bool> allowedVertex,
                                                  std::vector<double> primes,
                                                  std::vector<bool> neighbourhood);

std::vector<double> primeCountFixedVertex(const std::vector<std::vector<double>> &adjacencyMatrix,
                                          int vertex,
                                          unsigned long length,
                                          std::vector<int> subgraph,
                                          unsigned long neighboursNumber,
                                          std::vector<double> primes);

std::vector<std::vector<double>> subgraphAdjacencyMatrix(const std::vector<std::vector<double>> &adjacencyMatrix,
                                                         const std::vector<int> &subgraph);

double *subgraphAdjacencyArray(const std::vector<std::vector<double>> &adjacencyMatrix,
                               const std::vector<int> &subgraph);

std::vector<std::vector<double>> restrictedAdjacencyMatrix(const std::vector<std::vector<double>> &adjacencyMatrix,
                                                           const std::vector<int> &subgraph);

std::vector<std::vector<double>> matrixMultiplication(const double x,
                                                      const std::vector<std::vector<double>> &M);

std::vector<std::vector<double>> matrixAddition(const std::vector<std::vector<double>> &A,
                                                const std::vector<std::vector<double>> &B);

std::vector<std::vector<double>> squareMatrixMultiplication(const std::vector<std::vector<double>> &A,
                                                            const std::vector<std::vector<double>> &B);

long unsigned countTrue(const std::vector<bool> &vector);

double sum(double *array, long unsigned n);

double trace(const std::vector<std::vector<double>> &M);

#endif