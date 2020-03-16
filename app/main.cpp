#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "cycleCount.h"

int main() {
    std::vector<std::vector<double>> adjacencyMatrix = {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
                                                        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
                                                        {0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                                                        {0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                                                        {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0},
                                                        {0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0},
                                                        {0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0},
                                                        {0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0},
                                                        {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0},
                                                        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
                                                        {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0},
                                                        {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0},
                                                        {0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0},
                                                        {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                                                        {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                                                        {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}};

    // Compute simple cycles of length up to 10 in the graph
    const std::vector<double> countAll = cycleCount(adjacencyMatrix, 10, false);
    for (const auto c : countAll) {
        std::cout << std::setw(15) << std::left << round(c);
    }
    std::cout << std::endl;

    // Compute simple cycles of all length (length = adjacencyMatrix.size()) passing through vertex 3
    const std::vector<double> countFixedVertex = cycleCountFixedVertex(adjacencyMatrix, 3, adjacencyMatrix.size());
    for (const auto c : countFixedVertex) {
        std::cout << std::setw(15) << std::left << round(c);
    }
    std::cout << std::endl;
    return 0;
}

