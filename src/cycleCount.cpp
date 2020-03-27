#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "../header/eigenvalues.h"
#include "../header/cycleCount.h"

// EIGEN_MAX_ITERATION limit the max number of iteration to
// find eigen values to compute trace recursively
#define EIGEN_MAX_ITERATION 10000


/**
 * Counts all simple paths and cycles of length up to length included on adjacencyMatrix
 * @param adjacencyMatrix Undirected adjacency matrix of the graph G, this matrix may be weighted.
 * @param length maximum length of the simple path to be counted
 * @param directed true if the graph is directed, false if not
 * @return an array whose entry [i][x][y] is the number of simple paths of length i from x to yin the graph
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<std::vector<std::vector<double>>> pathsCount(std::vector<std::vector<double>> &adjacencyMatrix,
                                                         unsigned long length, bool directed) {

    // Gets rid of the self-loops
    std::vector<int> selfLoops(adjacencyMatrix.size(), 0);
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        selfLoops.at(i) = adjacencyMatrix.at(i).at(i);
        adjacencyMatrix.at(i).at(i) = 0;
    }
    unsigned long size = adjacencyMatrix.size(); // Number of vertices
    // Checks that maximum length does not exceed the total number of vertices
    length = length + 1;
    if (length > size) {
        length = size + 1;
    }
    // Initialisation of the simple path counting cell
    std::vector<double> line(size, 0);
    std::vector<std::vector<double>> row(size, line);
    std::vector<std::vector<std::vector<double>>> primes(length - 1, row);

    // Indicator vector of vertices that one may consider adding to a subgraph, initially all
    std::vector<bool> allowedVertex(size, true);

    // Loops over the vertices of the graph
    for (int vertexIndex = 0; vertexIndex < size; ++vertexIndex) {
        // Vertex visited is now forbidden
        allowedVertex.at(vertexIndex) = false;
        // Initialize indicator vector for neighbourhood
        std::vector<bool> neighbourhood(size, true);
        // current subgraph
        // vertices reachable via one edge
        for (int j = 0; j < neighbourhood.size(); ++j) {
            neighbourhood.at(j) = (adjacencyMatrix.at(vertexIndex).at(j) > 0);
        }
        neighbourhood.at(vertexIndex) = true;
        // Forms larger subgraphs containing the current one
        std::vector<int> subgraph = {vertexIndex};
        primes = recursiveSubgraphsPaths(adjacencyMatrix, length, subgraph,
                                         allowedVertex, primes, neighbourhood, directed);
    }

    // Puts self-loops back in the count
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        primes[0][i][i] += selfLoops[i];
        adjacencyMatrix[i][i] = selfLoops[i];
    }
    return primes;
}


/**
 * Finds all the connected induced subgraphs of size up "length" of a graph known through its adjacency matrix "A" and containing the subgraph "Subgraph"
 * @param adjacencyMatrix adjacency matrix of the graph, preferably sparse
 * @param length maximum subgraph size, an integer
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param allowedVertex indicator vector of pruned vertices that may be considered for addition to the current subgraph to form a larger one
 * @param primes list regrouping the contribution of all the subgraphs found so far
 * @param neighbourhood indicator vector of the vertices that are contained in the current subgraph or reachable via one edge
 * @param directed true if the graph is directed, false if not
 * @return primes with one more subgraph
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<std::vector<std::vector<double>>> recursiveSubgraphsPaths(std::vector<std::vector<double>> adjacencyMatrix,
                                                                      unsigned long length,
                                                                      std::vector<int> subgraph,
                                                                      std::vector<bool> allowedVertex,
                                                                      std::vector<std::vector<std::vector<double>>> primes,
                                                                      std::vector<bool> neighbourhood,
                                                                      bool directed) {

    // Counts the neighbours of the subgraph that are not contained in the subgraph itself
    long unsigned L = subgraph.size();
    long unsigned neighboursNumber = countTrue(neighbourhood) - L;

    // Gets the subgraph contribution to the prime list
    // Subgraphs of size 1 host no simple paths
    if (subgraph.size() > 1) {
        primes = primeCountPath(adjacencyMatrix, length, subgraph, neighboursNumber, primes);
    }
    if (L == length) {
        return primes;
    }
    // Indices of vertices that are both, in the neighbourhood and still allowed the vertices in the current subgraph are not allowed
    std::vector<int> neighbours;
    for (int i = 0; i < neighbourhood.size(); ++i) {
        if (neighbourhood[i] != 0 && allowedVertex.at(i) == true) {
            neighbours.push_back(i);
        }
    }

    // Adds each neighbour found above to Subgraph to form a new subgraph
    for (int v : neighbours) {
        // New subgraph
        if (subgraph.size() > L) {
            subgraph[L] = v;
        } else {
            subgraph.push_back(v);
        }
        // Vertex just added to new subgraph becomes forbidden
        allowedVertex.at(v) = false;
        // Calculate new neighbourhood
        std::vector<bool> newNeighbourhood(neighbourhood.size());
        for (int j = 0; j < neighbourhood.size(); ++j) {
            newNeighbourhood.at(j) = neighbourhood.at(j) || (adjacencyMatrix.at(v).at(j) > 0);
        }
        primes = recursiveSubgraphsPaths(adjacencyMatrix, length, subgraph, allowedVertex,
                                         primes, newNeighbourhood, directed);
    }
    return primes;
}


/**
 * Calculates the contribution to the combinatorial sieve of a given subgraph. This function is an implementation
 * of the equation extracting prime numbers from connected induced subgraphs.
 * @param adjacencyMatrix adjacency matrix of the graph, must be symmetric
 * @param length maximum subgraph size, an integer
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param neighboursNumber number of neighbours from the induced subgraph to the graph
 * @param primes current paths count array
 * @return updated paths count array
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<std::vector<std::vector<double>>> primeCountPath(const std::vector<std::vector<double>> &adjacencyMatrix,
                                                             unsigned long length,
                                                             std::vector<int> subgraph, unsigned long neighboursNumber,
                                                             std::vector<std::vector<std::vector<double>>> primes) {

    // Extracts the adjacency matrix x of the subgraph
    unsigned long size = adjacencyMatrix.size();
    unsigned long subgraphSize = subgraph.size();
    std::vector<std::vector<double>> x = subgraphAdjacencyMatrix(adjacencyMatrix, subgraph);
    std::vector<std::vector<double>> xPath = x;
    for (int i = 0; i < subgraphSize - 2; ++i) {
        xPath = squareMatrixMultiplication(xPath, x);
    }
    // Maximum value of k yielding a relevant non-zero Binomial(N(H),k-|H|)
    auto mk = std::min(length - 1, neighboursNumber + subgraphSize - 1);

    // The initial value of Binomial(N(H),k-|H|)
    unsigned long binomialCoeff = 1;
    for (unsigned long k = subgraphSize - 1; k <= mk; ++k) {
        // Update the simple path count at the right place (entry ij for paths from i to j)
        double pathValue = pow(-1.0, (double) k) * (double) binomialCoeff *
                           pow(-1.0, (double) subgraphSize - 1);
        std::vector<std::vector<double>> u = matrixMultiplication(pathValue, xPath);
        for (int i = 0; i < subgraphSize; ++i) {
            for (int j = 0; j < subgraphSize; ++j) {
                primes[k - 1][subgraph[i]][subgraph[j]] += u[i][j];
            }
        }

        // Update the simple cycle count at the right place (i.e. entry ii for cycles from i to itself)
        if (k + 1 <= length - 1) {
            std::vector<std::vector<double>> v1 = squareMatrixMultiplication(x, xPath);
            double cycleValue = pow(-1.0, (double) k + 1) * (double) binomialCoeff *
                                pow(-1.0, (double) subgraphSize);
            std::vector<std::vector<double>> v2 = matrixMultiplication(cycleValue, v1);
            for (int i = 0; i < subgraphSize; ++i) {
                primes[k][subgraph[i]][subgraph[i]] += v2[i][i];
            }
        }
        // Preparation of next loop
        xPath = squareMatrixMultiplication(xPath, x);
        // Increase by one value of k in a binomial coefficient
        binomialCoeff = binomialCoeff * (subgraphSize - (k + 1) + neighboursNumber) / (1 - subgraphSize + (k + 1));

    }
    return primes;
}


/**
 * Counts all simple cycle of length up to length included on adjacencyMatrix
 * @param adjacencyMatrix Undirected adjacency matrix of the graph G, this matrix may be weighted.
 * @param length maximum length of the simple cyles to be counted
 * @param directed true if the graph is directed, false if not
 * @return an array whose entry i is the number of simple cycles of length i in the graph
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> cycleCount(std::vector<std::vector<double>> &adjacencyMatrix,
                               unsigned long length, bool directed) {

    // Gets rid of the self-loops
    std::vector<int> selfLoops(adjacencyMatrix.size(), 0);
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        selfLoops.at(i) = adjacencyMatrix.at(i).at(i);
        adjacencyMatrix.at(i).at(i) = 0;
    }
    unsigned long size = adjacencyMatrix.size(); // Number of vertices
    // Checks that maximum length does not exceed the total number of vertices
    if (length > size) {
        length = size;
    }
    // Initialisation of the simple path counting cell
    std::vector<double> primes(length, 0);

    // Indicator vector of vertices that one may consider adding to a subgraph, initially all
    std::vector<bool> allowedVertex(size, true);

    // Loops over the vertices of the graph
    for (int vertexIndex = 0; vertexIndex < size; ++vertexIndex) {
        // Vertex visited is now forbidden
        allowedVertex.at(vertexIndex) = false;
        // Initialize indicator vector for neighbourhood
        std::vector<bool> neighbourhood(size, true);
        // current subgraph
        // vertices reachable via one edge
        for (int j = 0; j < neighbourhood.size(); ++j) {
            neighbourhood.at(j) = (adjacencyMatrix.at(vertexIndex).at(j) > 0);
        }
        neighbourhood.at(vertexIndex) = true;
        // Forms larger subgraphs containing the current one
        std::vector<int> subgraph = {vertexIndex};
        primes = recursiveSubgraphs(adjacencyMatrix, length, subgraph,
                                    allowedVertex, primes, neighbourhood, directed);
    }

    // Puts self-loops back in the count
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        primes[0] += selfLoops[i];
        adjacencyMatrix[i][i] = selfLoops[i];
    }
    return primes;
}

/**
 * Finds all the connected induced subgraphs of size up "length" of a graph known through its adjacency matrix "A" and containing the subgraph "Subgraph"
 * @param adjacencyMatrix adjacency matrix of the graph, preferably sparse
 * @param length maximum subgraph size, an integer
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param allowedVertex indicator vector of pruned vertices that may be considered for addition to the current subgraph to form a larger one
 * @param primes list regrouping the contribution of all the subgraphs found so far
 * @param neighbourhood indicator vector of the vertices that are contained in the current subgraph or reachable via one edge
 * @param directed true if the graph is directed, false if not
 * @return primes with one more subgraph
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> recursiveSubgraphs(std::vector<std::vector<double>> adjacencyMatrix,
                                       unsigned long length,
                                       std::vector<int> subgraph,
                                       std::vector<bool> allowedVertex,
                                       std::vector<double> primes,
                                       std::vector<bool> neighbourhood,
                                       bool directed) {

    // Counts the neighbours of the subgraph that are not contained in the subgraph itself
    long unsigned L = subgraph.size();
    long unsigned neighboursNumber = countTrue(neighbourhood) - L;
    // Gets the subgraph contribution to the prime list
    primes = primeCount(adjacencyMatrix, length, subgraph, neighboursNumber, primes, directed);
    if (L == length) {
        return primes;
    }
    // Indices of vertices that are both, in the neighbourhood and still allowed the vertices in the current subgraph are not allowed
    std::vector<int> neighbours;
    for (int i = 0; i < neighbourhood.size(); ++i) {
        if (neighbourhood[i] != 0 && allowedVertex.at(i) == true) {
            neighbours.push_back(i);
        }
    }

    // Adds each neighbour found above to Subgraph to form a new subgraph
    for (int v : neighbours) {
        // New subgraph
        if (subgraph.size() > L) {
            subgraph[L] = v;
        } else {
            subgraph.push_back(v);
        }
        // Vertex just added to new subgraph becomes forbidden
        allowedVertex.at(v) = false;
        // Calculate new neighbourhood
        std::vector<bool> newNeighbourhood(neighbourhood.size());
        for (int j = 0; j < neighbourhood.size(); ++j) {
            newNeighbourhood.at(j) = neighbourhood.at(j) || (adjacencyMatrix.at(v).at(j) > 0);
        }
        primes = recursiveSubgraphs(adjacencyMatrix, length, subgraph, allowedVertex,
                                    primes, newNeighbourhood, directed);
    }
    return primes;
}


/**
 * Simple handler to split between directed and undirected function calls.
 * @param adjacencyMatrix adjacency matrix of the graph, must be symmetric
 * @param length maximum subgraph size, an integer
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param neighboursNumber number of neighbours from the induced subgraph to the graph
 * @param primes current cycles count array
 * @param directed true if the graph is directed, false if not
 * @return updated cycles count array
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> primeCount(const std::vector<std::vector<double>> &adjacencyMatrix,
                               unsigned long length,
                               std::vector<int> subgraph,
                               unsigned long neighboursNumber,
                               std::vector<double> primes,
                               bool directed) {
    if (directed) {
        return primeCountUndirected(adjacencyMatrix, length, subgraph, neighboursNumber, primes);
    } else {
        return primeCountDirected(adjacencyMatrix, length, subgraph, neighboursNumber, primes);
    }
}


/**
 * Calculates the contribution to the combinatorial sieve of a given subgraph. This function is
 * an implementation of the equation extracting prime numbers from connected induced subgraphs.
 * It use eigenvalues to compute trace to the power recursively.
 * @protected do not use externally, use primeCount(_, false) instead
 * @param adjacencyMatrix adjacency matrix of the graph, must be symmetric
 * @param length maximum subgraph size, an integer
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param neighboursNumber number of neighbours from the induced subgraph to the graph
 * @param primes current cycles count array
 * @return updated cycles count array
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> primeCountUndirected(const std::vector<std::vector<double>> &adjacencyMatrix, unsigned long length,
                                         std::vector<int> subgraph, unsigned long neighboursNumber,
                                         std::vector<double> primes) {

    unsigned long subgraphSize = subgraph.size();
    // Extracts the adjacency matrix subgraphMatrix of the subgraph
    double *xDoubleArray = subgraphAdjacencyArray(adjacencyMatrix, subgraph);
    //  Get list of eigenvalues, each to power |H| to later compute Trace(A_H^k) recursively
    auto *eigenVectors = new double[subgraphSize * subgraphSize];
    auto eigenValues = new double[subgraphSize];
    int iterationsCount, rotationCount;
    computeEigenvalues(subgraphSize, xDoubleArray, EIGEN_MAX_ITERATION, eigenVectors,
                       eigenValues, iterationsCount, rotationCount);
    auto xS = new double[subgraphSize];
    for (int i = 0; i < subgraphSize; ++i) {
        xS[i] = pow(eigenValues[i], subgraphSize);
    }
    //Maximum value of k yielding a relevant non-zero Binomial(N(H),k-|H|)
    auto mk = std::min(length, neighboursNumber + subgraphSize);

    // The initial value of Binomial(N(H),k-|H|)
    unsigned long binomialCoeff = 1;
    for (unsigned long k = subgraphSize; k < mk; ++k) {
        // Combinatorial sieve
        primes[k - 1] += pow(-1.0, k) / (double) k * (double) binomialCoeff *
                         pow(-1.0, subgraphSize) * sum(xS, subgraphSize);

        // Preparation of next loop: list of eigenvalues to power k is put to power k+1 and value of k in a binomial coefficient is increased by 1
        for (int i = 0; i < subgraphSize; ++i) {
            xS[i] = eigenValues[i] * xS[i];
        }
        binomialCoeff = binomialCoeff * (subgraphSize - k + neighboursNumber) / (1 - subgraphSize + k);
    }

    // Contribution at maximum k value is separate to avoid incrementing xS, and the binomial coefficient one time too many
    primes[mk - 1] += pow(-1.0, mk) / (double) mk * (double) binomialCoeff *
                      pow(-1.0, subgraphSize) * sum(xS, subgraphSize);
    return primes;
}


/**
 * Calculates the contribution to the combinatorial sieve of a given subgraph. This function is an implementation
 * of the equation extracting prime numbers from connected induced subgraphs.
 * @protected do not use externally, use primeCount(_, true) instead
 * @param adjacencyMatrix adjacency matrix of the graph, must be symmetric
 * @param length maximum subgraph size, an integer
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param neighboursNumber number of neighbours from the induced subgraph to the graph
 * @param primes current cycles count array
 * @return updated cycles count array
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> primeCountDirected(const std::vector<std::vector<double>> &adjacencyMatrix, unsigned long length,
                                       std::vector<int> subgraph, unsigned long neighboursNumber,
                                       std::vector<double> primes) {

    // Extracts the adjacency matrix x of the subgraph
    unsigned long subgraphSize = subgraph.size();
    std::vector<std::vector<double>> initialX = subgraphAdjacencyMatrix(adjacencyMatrix, subgraph);
    // Compute Trace(A_H^k)
    std::vector<std::vector<double>> x = initialX;
    for (int i = 0; i < subgraphSize - 1; ++i) {
        x = squareMatrixMultiplication(x, initialX);
    }
    // Maximum value of k yielding a relevant non-zero Binomial(N(H),k-|H|)
    auto mk = std::min(length, neighboursNumber + subgraphSize);

    // The initial value of Binomial(N(H),k-|H|)
    unsigned long binomialCoeff = 1;
    for (unsigned long k = subgraphSize; k < mk; ++k) {
        // Combinatorial sieve
        primes[k - 1] += pow(-1.0, k) / (double) k * (double) binomialCoeff *
                         pow(-1.0, subgraphSize) * (double) trace(x);

        // Preparation of next loop: list of eigenvalues to power k is put to power k+1
        x = squareMatrixMultiplication(x, initialX);
        // Increase by one value of k in a binomial coefficient
        binomialCoeff = binomialCoeff * (subgraphSize - k + neighboursNumber) / (1 - subgraphSize + k);
    }

    // Contribution at maximum k value is separate to avoid incrementing xS, and the binomial coefficient one time too many
    primes[mk - 1] += pow(-1.0, mk) / (double) mk * (double) binomialCoeff *
                      pow(-1.0, subgraphSize) * (double) trace(x);
    return primes;
}


/**
 * Counts all simple cycle of length up to length passing through a specified vertex
 * @param adjacencyMatrix Undirected adjacency matrix of the graph G, this matrix may be weighted.
 * @param vertex the fixed vertex index
 * @param length maximum length of the simple cycles to be counted
 * @return an array whose entry i is the number of simple cycles of length i passing through the fixed vertex
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> cycleCountFixedVertex(std::vector<std::vector<double>> &adjacencyMatrix,
                                          int vertex,
                                          unsigned long length) {

    // Gets rid of the self-loops
    std::vector<int> selfLoops(adjacencyMatrix.size(), 0);
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        selfLoops.at(i) = adjacencyMatrix.at(i).at(i);
        adjacencyMatrix.at(i).at(i) = 0;
    }
    unsigned long size = adjacencyMatrix.size(); // Number of vertices
    // Checks that maximum length does not exceed the total number of vertices
    if (length > size) {
        length = size;
    }
    // Initialisation of the simple path counting cell
    std::vector<double> primes(length, 0);

    // Indicator vector of vertices that one may consider adding to a subgraph, initially all
    std::vector<bool> allowedVertex(length, true);

    // Loops over the vertices of the graph
    allowedVertex.at(vertex) = false;
    // Initialize indicator vector for neighbourhood
    std::vector<bool> neighbourhood(length, true);
    // Vertices reachable via one edge
    for (int j = 0; j < neighbourhood.size(); ++j) {
        neighbourhood.at(j) = (adjacencyMatrix.at(vertex).at(j) > 0);
    }
    neighbourhood.at(vertex) = true;
    // Forms larger subgraphs containing the current one
    std::vector<int> subgraph = {vertex};
    primes = recursiveSubgraphsFixedVertex(adjacencyMatrix, vertex, length, subgraph, allowedVertex, primes,
                                           neighbourhood);

    // Puts self-loops back in the count
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        primes[0] += selfLoops[i];
        adjacencyMatrix[i][i] = selfLoops[i];
    }
    return primes;
}

/**
 * Finds all the connected induced subgraphs of size up "length" of a graph known through its adjacency matrix "A" and containing the subgraph "Subgraph"
 * @param adjacencyMatrix adjacency matrix of the graph, must be symmetric
 * @param vertex the fixed vertex index
 * @param length maximum subgraph size, an integer
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param allowedVertex indicator vector of pruned vertices that may be considered for addition to the current subgraph to form a larger one
 * @param primes list regrouping the contribution of all the subgraphs found so far
 * @param neighbourhood indicator vector of the vertices that are contained in the current subgraph or reachable via one edge
 * @return updated cycles count list
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> recursiveSubgraphsFixedVertex(std::vector<std::vector<double>> adjacencyMatrix,
                                                  int vertex,
                                                  unsigned long length,
                                                  std::vector<int> subgraph,
                                                  std::vector<bool> allowedVertex,
                                                  std::vector<double> primes,
                                                  std::vector<bool> neighbourhood) {

    // Counts the neighbours of the subgraph that are not contained in the subgraph itself TODO : Bug here, check index
    long unsigned L = subgraph.size();
    long unsigned neighboursNumber = countTrue(neighbourhood) - L;
    // Gets the subgraph contribution to the prime list
    primes = primeCountFixedVertex(adjacencyMatrix, vertex, length, subgraph, neighboursNumber, primes);
    if (L == length) {
        return primes;
    }
    // Indices of vertices that are both, in the neighbourhood and still allowed the vertices in the current subgraph are not allowed
    std::vector<int> neighbours;
    for (int i = 0; i < neighbourhood.size(); ++i) {
        if (neighbourhood[i] != 0 && allowedVertex.at(i) == true) {
            neighbours.push_back(i);
        }
    }

    // Adds each neighbour found above to Subgraph to form a new subgraph
    for (int v : neighbours) {
        // New subgraph
        if (subgraph.size() > L) {
            subgraph[L] = v;
        } else {
            subgraph.push_back(v);
        }
        // Vertex just added to new subgraph becomes forbidden
        allowedVertex.at(v) = false;
        // Calculate new neighbourhood
        std::vector<bool> newNeighbourhood(neighbourhood.size());
        for (int j = 0; j < neighbourhood.size(); ++j) {
            newNeighbourhood.at(j) = neighbourhood.at(j) || (adjacencyMatrix.at(v).at(j) > 0);
        }
        primes = recursiveSubgraphsFixedVertex(adjacencyMatrix, vertex, length, subgraph, allowedVertex, primes,
                                               newNeighbourhood);
    }
    return primes;
}

/**
 * Calculates the contribution to the combinatorial sieve of a given subgraph. This function is an implementation of
 * the equation extracting prime numbers from connected induced subgraphs for a fixed vertex.
 * @param adjacencyMatrix adjacency matrix of the graph
 * @param length maximum subgraph size, an integer
 * @param vertex the fixed vertex
 * @param subgraph current subgraph, a list of vertices, further vertices are added to this list
 * @param neighboursNumber number of neighbours from the induced subgraph to the graph
 * @param primes current cycles count array
 * @return updated cycles count array
 * @remark Originally designed by P.-L. Giscard, N. Kriege, R. C. Wilson, July 2017
 */
std::vector<double> primeCountFixedVertex(const std::vector<std::vector<double>> &adjacencyMatrix,
                                          int vertex,
                                          unsigned long length, std::vector<int> subgraph,
                                          unsigned long neighboursNumber,
                                          std::vector<double> primes) {

    // Extracts the adjacency matrix x of the subgraph
    unsigned long subgraphSize = subgraph.size();
    std::vector<std::vector<double>> initialX = restrictedAdjacencyMatrix(adjacencyMatrix, subgraph);

    std::vector<std::vector<double>> x = initialX;
    for (int i = 0; i < subgraphSize - 1; ++i) {
        x = squareMatrixMultiplication(x, initialX);
    }
    // Maximum value of k yielding a relevant non-zero Binomial(N(H),k-|H|)
    auto mk = std::min(length, neighboursNumber + subgraphSize);

    // The initial value of Binomial(N(H),k-|H|)
    unsigned long binomialCoeff = 1;
    for (unsigned long k = subgraphSize; k < mk; ++k) {
        // Combinatorial sieve
        primes[k - 1] +=
                pow(-1.0, k) * (double) binomialCoeff * pow(-1.0, subgraphSize) * (double) x.at(vertex).at(vertex);

        // Preparation of next loop
        x = squareMatrixMultiplication(x, initialX);
        binomialCoeff = binomialCoeff * (subgraphSize - k + neighboursNumber) / (1 - subgraphSize + k);
    }

    // Contribution at maximum k value is separate to avoid incrementing xS, and the binomial coefficient one time too many
    primes[mk - 1] +=
            pow(-1.0, mk) * (double) binomialCoeff * pow(-1.0, subgraphSize) * (double) x.at(vertex).at(vertex);
    return primes;
}


/**
 * Get the adjacency matrix induced by a vertex list
 * @param adjacencyMatrix original adjacency matrix
 * @param subgraph list of index for the induced adjacency matrix
 * @return the adjacency matrix with only specified indexes
 */
std::vector<std::vector<double>> subgraphAdjacencyMatrix(const std::vector<std::vector<double>> &adjacencyMatrix,
                                                         const std::vector<int> &subgraph) {

    int n = subgraph.size();
    std::vector<std::vector<double>> x;
    x.reserve(n);
    for (int const &v1: subgraph) {
        std::vector<double> line;
        line.reserve(n);
        for (int const &v2: subgraph) {
            line.push_back(adjacencyMatrix[v1][v2]);
        }
        x.push_back(line);
    }
    return x;
}


/**
 * Get a flatten adjacency matrix induced by a vertex list
 * @param adjacencyMatrix original adjacency matrix
 * @param subgraph list of index for the induced adjacency matrix
 * @return the adjacency matrix with only specified indexes
 */
double *subgraphAdjacencyArray(const std::vector<std::vector<double>> &adjacencyMatrix,
                               const std::vector<int> &subgraph) {

    int n = subgraph.size();
    auto *array = new double[n * n];
    int i = 0;
    for (int const &v1: subgraph) {
        for (int const &v2: subgraph) {
            array[i] = adjacencyMatrix[v1][v2];
            i++;
        }
    }
    return array;
}

/**
 * The adjacency matrix restricted to subgraph
 * @param adjacencyMatrix the original adjacency matrix
 * @param subgraph list of vertexes
 * @return adjacencyMatrix with (entry[i][j] = 0) if i or j not in subgraph
 */
std::vector<std::vector<double>> restrictedAdjacencyMatrix(const std::vector<std::vector<double>> &adjacencyMatrix,
                                                           const std::vector<int> &subgraph) {

    std::vector<std::vector<double>> x;
    x.reserve(adjacencyMatrix.size());
    for (int i = 0; i < adjacencyMatrix.size(); ++i) {
        std::vector<double> line(adjacencyMatrix.size(), 0);
        x.push_back(line);
    }
    for (int const &v1: subgraph) {
        for (int const &v2: subgraph) {
            x.at(v1).at(v2) = adjacencyMatrix.at(v1).at(v2);
        }
    }
    return x;
}

/**
 * @param x double value for multiplication
 * @param M square matrix
 * @return The matrix (x * M)
 */
std::vector<std::vector<double>> matrixMultiplication(const double x,
                                                      const std::vector<std::vector<double>> &M) {
    int n1 = M.size();
    int n2 = M[0].size();
    std::vector<std::vector<double>> MM(n1, std::vector<double>(n2, 0));
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            MM[i][j] += x * M[i][j];
        }
    }
    return MM;
}


/**
 * @param A first square matrix
 * @param B second square matrix
 * @return The matrix (A + B)
 */
std::vector<std::vector<double>> matrixAddition(const std::vector<std::vector<double>> &A,
                                                const std::vector<std::vector<double>> &B) {
    int n1 = A.size();
    int n2 = A[0].size();
    std::vector<std::vector<double>> M(n1, std::vector<double>(n2, 0));
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            M[i][j] += A[i][j] + B[i][j];
        }
    }
    return M;
}

/**
 * @param A first square matrix
 * @param B second square matrix
 * @return The matrix (A * B)
 */
std::vector<std::vector<double>> squareMatrixMultiplication(const std::vector<std::vector<double>> &A,
                                                            const std::vector<std::vector<double>> &B) {

    int n = A.size();
    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0));
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            for (int i = 0; i < n; ++i) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

/**
 * Count number of true entry in a vector
 * @param vector the boolean vector
 * @return true values counter
 */
long unsigned countTrue(const std::vector<bool> &vector) {
    long unsigned count = 0;
    for (int i : vector) {
        if (i) {
            count++;
        }
    }
    return count;
}

/**
 * @param array an array of double
 * @param n an index
 * @return Sum of the n first values of array
 */
double sum(double *array, long unsigned n) {
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += array[i];
    }
    return sum;
}

/**
 * @param M an array of double
 * @return Trace of the matrix
 */
double trace(const std::vector<std::vector<double>> &M) {
    double trace = 0;
    for (int i = 0; i < M.size(); ++i) {
        trace += M.at(i).at(i);
    }
    return trace;
}
