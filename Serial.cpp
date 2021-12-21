#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <map>
#include <set>
#include <chrono>

int INF = std::numeric_limits<int>::max();

struct Node {
    Node(int k, int c) {
        key = k;
        cost = c;
    }

    int key;
    int cost;
};

struct compare {
public:
    bool operator()(Node *a, Node *b) {

        return a->cost > b->cost;
    }
};

std::map<int, std::vector<Node *>> createAdjList(const std::string &filename) {
    std::fstream file(filename);
    int vertices;
    file >> vertices;

    std::map<int, std::vector<Node *>> adjList;
    for (int i = 0; i < vertices; ++i) {
        std::vector<Node *> tempVec;
        for (int j = 0; j < vertices; ++j) {
            std::string temp;
            file >> temp;
            int cost = std::stoi(temp);
            if (cost != 0) {
                tempVec.push_back(new Node(j, cost));
            }
        }
        adjList[i] = tempVec;
    }
    return adjList;
}

int getSourceNodeFromAdjList(std::map<int, std::vector<Node *>> &adjList) {
    for (const auto &i: adjList) {
        if (!i.second.empty()) {
            return i.first;
        }
    }
    return -1;
}

void printAdjList(const std::map<int, std::vector<Node *>> &adjList) {
    for (const auto &i: adjList) {
        std::cout << i.first << ":\t";
        for (auto j: i.second) {
            std::cout << j->key << " -> " << j->cost << "\t";
        }
        std::cout << std::endl;
    }
}

void DijkstraAdjList(std::map<int, std::vector<Node *>> &adjList, int source) {
    std::vector<int> dist(adjList.size(), INF);
    std::priority_queue<Node *, std::vector<Node *>, compare> PQ;
    std::set<int> PQList;
    dist[source] = 0;
    PQ.push(new Node(source, dist[source]));
    PQList.insert(source);

    while (!PQ.empty()) {
        auto u = PQ.top();
        PQ.pop();
        PQList.erase(u->key);
        auto neighbours = adjList[u->key];
        for (auto n: neighbours) {
            bool is_in_PQ = PQList.find(n->key) != PQList.end();
            int alt = dist[u->key] + n->cost;
            if (alt < dist[n->key]) {
                dist[n->key] = alt;
                if (!is_in_PQ) {
                    PQ.push(n);
                    PQList.insert(n->key);
                }
            }
        }
    }
    std::cout << "Distances" << std::endl;
    for (int i = 0; i < adjList.size(); ++i) {
        if (i != source) {
            std::cout << i << ": " << dist[i] << std::endl;
        }

    }
}

std::map<int, std::vector<Node *>> MatrixToList(std::vector<std::vector<int>> adjMatrix) {
    std::map<int, std::vector<Node *>> adjList;
    for (int i = 0; i < adjMatrix.size(); ++i) {
        std::vector<Node *> temp;
        for (int j = 0; j < adjMatrix[i].size(); ++j) {
            int cost = adjMatrix[i][j];
            if (cost != 0) {
                temp.push_back(new Node(j, adjMatrix[i][j]));
            }
        }
        adjList[i] = temp;
    }
    return adjList;
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

void DijkstraPQ(std::vector<std::vector<int>> adjMatrix, int src) {
    std::vector<int> dist(adjMatrix.size(), INF);
    std::priority_queue<std::pair<int, int>> PQ;
    std::set<int> PQList;
    dist[src] = 0;
    PQ.push(std::make_pair(src, dist[src]));
    PQList.insert(src);
    while (!PQ.empty()) {
        auto u = PQ.top();
        PQ.pop();
        PQList.erase(u.first);
        auto neighbours = adjMatrix[u.first];
        for (int i = 0; i < neighbours.size(); ++i) {
            //i is the key -> neighbour[i] is the cost
            if (neighbours[i] == 0) {
                continue;
            }
            bool is_in_PQ = PQList.find(i) != PQList.end();
            int alt = dist[u.first] + neighbours[i];
            if (alt < dist[i]) {
                dist[i] = alt;
                if (!is_in_PQ) {
                    PQ.push(std::make_pair(i, neighbours[i]));
                    PQList.insert(i);
                }
            }
        }
    }
    std::cout << "Distances" << std::endl;
    for (int i = 0; i < adjMatrix.size(); ++i) {
        if (i != src) {
            std::cout << i << ": " << dist[i] << std::endl;
        }

    }
}

std::vector<int> Dijkstra(std::vector<std::vector<int>> adjMatrix, int src) {
    int N = (int) adjMatrix.size();
    int count = 0;

    std::vector<int> dist(N, INF);
    std::vector<bool> visited(N, false);

    for (int v = 0; v < N; v++) {
        if (v == src) {
            continue;
        }
        int cost = adjMatrix[src][v];
        if (cost != 0) {
            dist[v] = cost;
        }

    }
    visited[src] = true;
    count++;


    for (int i = 0; i < N; i++) {
        int minDist = INF;
        int u;
        for (int v = 0; v < N; v++) {
            if (!visited[v]) {
                if (dist[v] < minDist) {
                    minDist = dist[v];
                    u = v;
                }
            }
        }

        visited[u] = true;
        count++;

        for (int v = 0; v < N; v++) {
            if (!visited[v]) {
                if (adjMatrix[u][v] != 0) {
                    int newDist = minDist + adjMatrix[u][v];
                    dist[v] = std::min(dist[v], newDist);
                }

            }
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
            output << i << ": " << "Isolated"<<std::endl;
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
    std::string outputName = testName + "_Serial_Output.txt";

    auto adjMatrix = createAdjMatrix(fileName);

    int src = getSourceNodeFromAdjMatrix(adjMatrix);
    auto start = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch());
    auto dist = Dijkstra(adjMatrix, src);

    auto end = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch());
    auto diff = end - start;
    writeToFile(dist, src, outputName);

    std::ofstream results("Results.txt", std::ios_base::app);
    results << "Serial Time Taken:\t\t\t" << diff.count() << std::endl;
    std::cout << outputName;


    return 0;
}