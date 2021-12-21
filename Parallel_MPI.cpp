#include <iostream>
#include <mpi.h>
#include <string>
#include <fstream>
#include <vector>
#include <thread>
#include <cmath>

int INF = std::numeric_limits<int>::max();

int getSourceNode(const int adjMatrix[], int size) {
    int N = (int) sqrt(size);
    int idx = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N && idx < N * N; ++j) {
            if (adjMatrix[idx] == 0) {
                return i;
            }
            idx++;
        }
    }
    return -1;
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

MPI_Datatype createDataType(int vertices, int localVertices) {
    MPI_Aint lb, extent;
    MPI_Datatype block;
    MPI_Datatype firstBlockCol;
    MPI_Datatype subCol;

    MPI_Type_contiguous(localVertices, MPI_INT, &block);
    MPI_Type_get_extent(block, &lb, &extent);

    MPI_Type_vector(vertices, localVertices, vertices, MPI_INT, &firstBlockCol);

    MPI_Type_create_resized(firstBlockCol, lb, extent, &subCol);

    MPI_Type_commit(&subCol);

    MPI_Type_free(&block);
    MPI_Type_free(&firstBlockCol);
    return subCol;
}


int main(int argc, char *argv[]) {
    double myTime;
    int size, rank;
    int vertices, source, localVertices;
    int info[3];
    int *adjMatrix;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fstream file;

    if (rank == 0) {
        if (argc < 2) {
            std::cout << "Missing Input File Argument" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return -1;
        }


        std::string fileName = argv[1];
        file.open(fileName);

        file >> vertices;
        if (vertices % size != 0) {
            std::cout << "Vertices Must Be Divisible By Number of Processors " << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return -1;
        }

        localVertices = vertices / size;
        int gridSize = vertices * vertices;
        adjMatrix = new int[gridSize];
        for (int i = 0; i < gridSize; ++i) {
            std::string temp;
            file >> temp;
            adjMatrix[i] = (std::stoi(temp));
        }
        source = getSourceNode(adjMatrix, gridSize);
        info[0] = vertices;
        info[1] = localVertices;
        info[2] = source;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    myTime = MPI_Wtime();



    MPI_Bcast(info, 3, MPI_INT, 0, MPI_COMM_WORLD);
    vertices = info[0];
    localVertices = info[1];
    source = info[2];

    int Dist[vertices];
    int localDist[localVertices];
    int localAdjMatrix[localVertices * vertices];
    int firstVertex = rank * localVertices;
    int lastVertex = firstVertex + localVertices - 1;

    MPI_Datatype subCol = createDataType(vertices, localVertices);
    MPI_Scatter(adjMatrix, 1, subCol, localAdjMatrix, localVertices * vertices, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<bool> marked(localVertices, true);

    int localMin[2];
    int globalMin[2];

    for (int j = 0; j < localVertices; j++) {
        int cost = localAdjMatrix[source * localVertices + j];
        if (cost != 0) {
            localDist[j] = cost;
        } else {
            localDist[j] = INF;
        }
    }

    for (int i = 0; i < localVertices; ++i) {
        marked[i] = true;
    }

    if (source >= firstVertex && source <= lastVertex) {
        marked[source - firstVertex] = false;
        localDist[source - firstVertex] = 0;
    }

    for (int i = 0; i < vertices; ++i) {
        localMin[0] = INF;
        localMin[1] = -1;
        for (int j = 0; j < localVertices; ++j) {
            if (marked[j] && localDist[j] < localMin[0]) {
                localMin[0] = localDist[j];
                localMin[1] = firstVertex + j;
            }
        }
        MPI_Allreduce(localMin, globalMin, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        int minDist = globalMin[0];
        int u = globalMin[1];

        if (u == localMin[1]) {
            marked[u - firstVertex] = false;
        }
        for (int j = 0; j < localVertices; j++) {
            if (marked[j]) {
                if (localAdjMatrix[u * localVertices + j] != 0) {
                    int newDist = minDist + localAdjMatrix[u * localVertices + j];
                    localDist[j] = std::min(localDist[j], newDist);
                }

            }
        }
    }

    MPI_Gather(localDist, localVertices, MPI_INT, Dist, localVertices, MPI_INT, 0, MPI_COMM_WORLD);


    myTime = MPI_Wtime() - myTime;
    double maxTime;
    MPI_Reduce(&myTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        std::string fullname = argv[1];
        std::string testName = fullname.substr(0, fullname.find_last_of('.'));
        testName = testName.substr(testName.find_last_of("/\\") + 1);
        std::string outputName = testName + "_MPI_Output.txt";

        std::ofstream output(outputName);
        output << "Distances" << std::endl;
        for (int i = 0; i < vertices; ++i) {
            if (i == source) {
                continue;
            }
            if (Dist[i] == INF || Dist[i] == -INF) {
                output << i << ": " << "Isolated" << std::endl;
                continue;
            }
            output << i << ": " << Dist[i] << std::endl;
        }
        output.close();

        std::ofstream results("Results.txt", std::ios_base::app);
        results << "MPI Parallel Time Taken:\t" << maxTime << std::endl;
        std::cout << outputName;
    }
    MPI_Finalize();
    return 0;
}