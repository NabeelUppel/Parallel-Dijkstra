#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <omp.h>
#include <sstream>
#include <limits>

int INF = std::numeric_limits<int>::max();

void printDist(std::vector<int> dist) {
    for (int i = 0; i < dist.size(); ++i) {
        std::cout << i << ": " << dist[i] << std::endl;
    }
    std::cout << std::endl << std::endl;
}

void printAdjMatrix(const std::vector<std::vector<int>> &adjMatrix) {
    for (auto &i: adjMatrix) {
        for (int j: i) {
            std::cout << j << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
}

std::vector<int> toIntRow(const std::string &line) {
    std::string tmp;
    std::vector<int> row;
    std::stringstream ss(line);

    while (getline(ss, tmp, ' ')) {
        row.push_back(std::stoi(tmp));
    }
    return row;
}


std::vector<std::vector<int>> createAdjMatrix(const std::string &filename) {
    std::fstream file(filename);
    int vertices;
    file >> vertices;
    std::vector<std::vector<int>> adjMatrix(vertices, std::vector<int>(vertices));
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            std::string temp;
            file >> temp;
            int cost = std::stoi(temp);
            adjMatrix[i][j] = cost;
        }
    }
    return adjMatrix;
}


int getSourceNodeFromAdjMatrix(std::vector<std::vector<int>> adjMatrix) {
    for (int i = 0; i < adjMatrix.size(); ++i) {
        for (int j = 0; j < adjMatrix[i].size(); ++j) {
            if (adjMatrix[i][j] != 0) {
                return i;
            }
        }
    }
    return -1;
}

std::vector<int> Dijkstra(std::vector<std::vector<int>> adjMatrix, int src) {
    int N = (int) adjMatrix.size();
    std::vector<int> dist(N, INF);
    std::vector<bool> visited(N, false);

    for (int v = 0; v < N; v++) {
        int cost = adjMatrix[src][v];
        if (cost != 0) {
            dist[v] = cost;
        }

    }

    visited[src] = true;
    dist[src] = 0;
    int minDist, u, localMinDist, localU, id;
    int firstVertex, lastVertex, thread_count;
#pragma omp parallel  private(id, firstVertex, lastVertex, localMinDist, localU, thread_count) shared (visited, minDist, dist, u, adjMatrix, N)
    {
        id = omp_get_thread_num();
        thread_count = omp_get_num_threads();
        firstVertex = (id * N) / thread_count;
        lastVertex = ((id + 1) * N) / thread_count - 1;
        for (int i = 0; i < N; i++) {
#pragma omp single
            {
                minDist = INF;
                u = -1;
            }
            localMinDist = INF;
            localU = -1;
            for (int k = firstVertex; k <= lastVertex; k++) {
                if (!visited[k]) {
                    if (dist[k] < localMinDist) {
                        localMinDist = dist[k];
                        localU = k;
                    }
                }
            }

#pragma omp critical
            {
                if (localMinDist < minDist) {
                    minDist = localMinDist;
                    u = localU;
                }

            }

# pragma omp barrier


# pragma omp single

            {
                if (u != -1) {
                    visited[u] = true;
                }
            }
# pragma omp barrier

            if (u != -1) {
                for (int v = firstVertex; v <= lastVertex; v++) {
                    if (!visited[v]) {
                        if (adjMatrix[u][v] != 0) {
                            int newDist = minDist + adjMatrix[u][v];
                            dist[v] = std::min(dist[v], newDist);
                        }
                    }

                }
            }

# pragma omp barrier
        }
    }
    return dist;
}

void writeToFile(std::vector<int> dist, int src, const std::string &fileName) {
    std::ofstream output(fileName);
    output << "Distances" << std::endl;
    for (int i = 0; i < dist.size(); ++i) {
        if (i == src) {
            continue;
        }
        if (dist[i] == INF || dist[i] == -INF) {
            output << i << ": " << "Isolated" << std::endl;
            continue;
        }
        output << i << ": " << dist[i] << std::endl;

    }
    output.close();
}


int main(int argc, char *argv[]) {
     if (argc < 2) {
        std::cout << "Missing Input File Argument" << std::endl;
        return -1;
    }
    std::string fileName = argv[1];
    std::string testName = fileName.substr(0, fileName.find_last_of('.'));
    testName = testName.substr(testName.find_last_of("/\\") + 1);
    std::string outputName = testName + "_OMP_Output.txt";
    auto adjMatrix = createAdjMatrix(fileName);
    int src = getSourceNodeFromAdjMatrix(adjMatrix);

    auto start = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch());
    auto dist = Dijkstra(adjMatrix, src);
    auto end = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch());
    auto diff = end - start;

    writeToFile(dist, src, outputName);
    std::ofstream results("Results.txt", std::ios_base::app);
    results << "OMP Parallel Time Taken:\t" << diff.count() << std::endl;
    std::cout << outputName;


    return 0;
}