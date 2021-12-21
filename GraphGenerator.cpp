#include <iostream>
#include <vector>
#include <cmath>
#include<fstream>
#include <random>

int calculateEdges(int vertices, double target) {
    int edges = (int) std::ceil((vertices * (vertices - 1)) * target) / 2;
    return edges;
}


void saveMatrix(const std::vector<std::vector<int>> &graph, const std::string &filename) {
    std::ofstream file(filename);
    file << graph.size() << std::endl;
    for (auto &i: graph) {
        for (int j: i) {
            file << j << " ";
        }
        file << std::endl;
    }
    file.close();
}


std::vector<std::vector<int>> generateGraph(int vertices) {
    double targetDensity = 0.35;
    int edges = calculateEdges(vertices, targetDensity);
    std::vector<int> rows(vertices, 4);
    std::vector<std::vector<int>> adjMatrix(vertices, std::vector<int>(vertices));
    int totalEdges = 0;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> weights(1, 25); // define the range
    std::uniform_int_distribution<> randIndex(0, vertices - 1); // define the range

    while (totalEdges < edges) {
        int i = randIndex(gen);
        int j = randIndex(gen);
        int randWeight = weights(gen);
        if (i != j && adjMatrix[i][j] == 0) {
            adjMatrix[i][j] = randWeight;
            adjMatrix[j][i] = randWeight;
            totalEdges++;
        }
    }
    return adjMatrix;
}

int main(int argc, char *argv[]) {
    int vertices;
    if (argc < 2) {
        std::cout << "./Executable <number of vertices>";
        return -1;
    }
    vertices = std::stoi(argv[1]);
    if (vertices < 3) {
        std::cout << "Number of Vertices is Too Small";
        return -1;
    }
    auto adjMatrix = generateGraph(vertices);
    std::string output = argv[1];
    output += +"_Vertex_Graph.txt";

    saveMatrix(adjMatrix, output);
    std::cout<<output;
    return 0;
}
