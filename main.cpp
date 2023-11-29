#include <iostream>
#include <vector>
#include <fstream>
#include <queue>

auto dejkstra(std::vector<std::vector<int64_t>> &gr1) {
    int n = gr1.size();
    std::vector<std::vector<std::pair<int, int64_t>>> gr(n);
    int a, b, c;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (gr1[i][j] != 0) {
                gr[i].push_back({j, gr1[i][j]});
            }
        }
    }
    auto start = std::chrono::high_resolution_clock::now();
    int s = 0;
    int64_t inf = std::numeric_limits<int64_t>::max();
    std::vector<int64_t> dist(n, inf);
    dist[s] = 0;
    std::vector<char> used(n);
    for (int i = 0; i < n; ++i) {
        int v = -1;
        for (int j = 0; j < n; ++j) {
            if (!used[j] && (v == -1 || dist[j] < dist[v])) {
                v = j;
            }
        }
        if (dist[v] == inf) {
            break;
        }
        used[v] = true;

        for (const auto &edge: gr[v]) {
            int u = edge.first;
            int64_t weight = edge.second;
            if (dist[v] + weight < dist[u]) {
                dist[u] = dist[v] + weight;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now() - start;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end).count();
}

auto johnson(std::vector<std::vector<std::pair<int, int64_t>>>& gr) {
    auto start = std::chrono::high_resolution_clock::now();
    int n = gr.size();
    std::vector<std::vector<int64_t>> dist(n, std::vector<int64_t>(n));
    std::vector<int64_t> h(n, 0);
    gr.resize(n + 1);
    for (int i = 0; i < n; ++i) {
       
        gr[n].emplace_back(i, 0);
    }

    for (int k = 0; k < n; ++k) {
        std::vector<int64_t> dist_bellman(n + 1, std::numeric_limits<int64_t>::max());
        dist_bellman[n] = 0;
        for (int i = 0; i < n; ++i) {
            for (const auto& edge : gr[i]) {
                int v = edge.first;
                int64_t weight = edge.second;
                dist_bellman[v] = std::min(dist_bellman[v], dist_bellman[i] + weight);
            }
        }
        h = dist_bellman;
    }

    for (int s = 0; s < n; ++s) {
        std::priority_queue<std::pair<int64_t, int>, std::vector<std::pair<int64_t, int>>, std::greater<std::pair<int64_t, int>>> pq;
        std::vector<int64_t> dist_dijkstra(n, std::numeric_limits<int64_t>::max());
        dist_dijkstra[s] = 0;
        std::vector<char> used(n);

        pq.emplace(0, s);

        while (!pq.empty()) {
            int v = pq.top().second;
            pq.pop();

            if (used[v]) {
                continue;
            }

            used[v] = true;

            for (const auto& edge : gr[v]) {
                int u = edge.first;
                int64_t weight = edge.second + h[v] - h[u];
                if (dist_dijkstra[v] + weight < dist_dijkstra[u]) {
                    dist_dijkstra[u] = dist_dijkstra[v] + weight;
                    pq.emplace(dist_dijkstra[u], u);
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            if (dist_dijkstra[i] != std::numeric_limits<int64_t>::max()) {
                dist[s][i] = dist_dijkstra[i] - h[s] + h[i];
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now() - start;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end).count();
}

auto floyduorshil(std::vector<std::vector<int64_t>> &dist1, int size) {
    std::vector<std::vector<int64_t>> dist(size);
    std::copy(dist1.begin(), dist1.begin() + size, dist.begin());
    auto start = std::chrono::high_resolution_clock::now();
    int n = dist.size();
    int64_t inf = std::numeric_limits<int64_t>::max();
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist[i][k] != inf && dist[k][j] != inf &&
                    dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now() - start;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end).count();
}

auto fordbellman(std::vector<std::pair<int64_t, std::pair<int, int>>> &gr, int size) {
    auto start = std::chrono::high_resolution_clock::now();
    int s = 0;
    int64_t inf = 21474836470;
    std::vector<int64_t> dist1(size, inf);
    dist1[s] = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < gr.size(); j++) {
            if (dist1[gr[j].second.first] < inf) {
                if (dist1[gr[j].second.second] > dist1[gr[j].second.first] + gr[j].first) {
                    dist1[gr[j].second.second] = dist1[gr[j].second.first] + gr[j].first;
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now() - start;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end).count();
}

std::vector<std::vector<int64_t>> generateFullGraph(int numVertices) {
    srand(time(nullptr));
    std::vector<std::vector<int64_t>> graph(numVertices, std::vector<int64_t>(numVertices, 0));
    for (int i = 0; i < numVertices; ++i) {
        for (int j = i + 1; j < numVertices; ++j) {
            int64_t weight = rand() % 10 + 1;
            graph[i][j] = weight;
            graph[j][i] = weight;
        }
    }
    return graph;
}

std::vector<std::vector<int64_t>> generateSecondGraph(int numVertices) {
    srand(time(nullptr));
    std::vector<std::vector<int64_t>> graph(numVertices, std::vector<int64_t>(numVertices, 0));
    for (int i = 1; i < numVertices; ++i) {
        int64_t weight = rand() % 10 + 1;
        graph[i][i - 1] = weight;
        graph[i - 1][i] = weight;
    }
    int extra_edges = numVertices / 2;
    for (int i = 0; i < extra_edges; ++i) {
        int u = rand() % numVertices;
        int v = rand() % numVertices;
        if (u != v) {
            int weight = rand() % 10 + 1;
            graph[u][v] = weight;
            graph[v][u] = weight;
        }
    }
    return graph;
}

std::vector<std::vector<int64_t>> generateThirdGraph(int numVertices) {
    srand(time(nullptr));
    std::vector<std::vector<int64_t>> graph(numVertices, std::vector<int64_t>(numVertices, 0));

    for (int i = 1; i < numVertices; ++i) {
        int64_t weight = rand() % 10 + 1;
        graph[i - 1][i] = weight;
        graph[i][i - 1] = weight;
    }

    return graph;
}

int cntEdg(std::vector<std::vector<int64_t>> &gr, int size) {
    int cnt = 0;
    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if (gr[i][j] != 0) {
                cnt++;
            }
        }
    }
    return cnt;
}

void start_dejkstra() {
    std::ofstream myfile;
    myfile.open("dejkstra.csv");
    myfile
            << "; full graph; second type graph; third type graph; full graph edg; second type graph edg; third type graph edg; \n";


    for (int i = 10; i <= 1010; i += 50) {
        std::vector<std::vector<int64_t>> fullGraph = generateFullGraph(i);
        std::vector<std::vector<int64_t>> secondType = generateSecondGraph(i);
        std::vector<std::vector<int64_t>> thirdType = generateThirdGraph(i);
        myfile << i << "; ";
        long long ans1 = 0;
        long long ans2 = 0;
        long long ans3 = 0;
        for (int j = 0; j < 5; j++) {
            ans1 += dejkstra(fullGraph);
            ans2 += dejkstra(secondType);
            ans3 += dejkstra(thirdType);
        }
        ans1 /= 5;
        ans2 /= 5;
        ans3 /= 5;
        myfile << ans1 << "; ";
        myfile << ans2 << "; ";
        myfile << ans3 << "; ";
        myfile << cntEdg(fullGraph, i) << "; ";
        myfile << cntEdg(secondType, i) << "; ";
        myfile << cntEdg(thirdType, i) << "; ";
        myfile << "\n";
    }
}

void start_johnson() {
    std::ofstream myfile;
    myfile.open("johnson.csv");
    myfile
            << "; full graph; second type graph; third type graph; full graph edg; second type graph edg; third type graph edg; \n";


    for (int i = 10; i <= 1010; i += 50) {
        std::vector<std::vector<int64_t>> fullGraph = generateFullGraph(i);
        std::vector<std::vector<int64_t>> secondType = generateSecondGraph(i);
        std::vector<std::vector<int64_t>> thirdType = generateThirdGraph(i);
        std::vector<std::vector<std::pair<int, int64_t>>> fullGraph1(i);
        std::vector<std::vector<std::pair<int, int64_t>>> secondType1(i);
        std::vector<std::vector<std::pair<int, int64_t>>> thirdType1(i);
        for (int i1 = 0; i1 < i; i1++) {
            for (int j = 0; j < i; j++) {
                if (fullGraph[i1][j] != 0) {
                    fullGraph1[i1].push_back({j, fullGraph[i1][j]});
                }
            }
        }
        for (int i1 = 0; i1 < i; i1++) {
            for (int j = 0; j < i; j++) {
                if (thirdType[i1][j] != 0) {
                    thirdType1[i1].push_back({j, thirdType[i1][j]});
                }
            }
        }
        for (int i1 = 0; i1 < i; i1++) {
            for (int j = 0; j < i; j++) {
                if (secondType[i1][j] != 0) {
                    secondType1[i1].push_back({j, secondType[i1][j]});
                }
            }
        }
        myfile << i << "; ";
        long long ans1 = 0;
        long long ans2 = 0;
        long long ans3 = 0;
        for (int j = 0; j < 5; j++) {
            ans1 += johnson(fullGraph1);
            ans2 += johnson(secondType1);
            ans3 += johnson(thirdType1);
        }
        ans1 /= 5;
        ans2 /= 5;
        ans3 /= 5;
        myfile << ans1 << "; ";
        myfile << ans2 << "; ";
        myfile << ans3 << "; ";
        myfile << cntEdg(fullGraph, i) << "; ";
        myfile << cntEdg(secondType, i) << "; ";
        myfile << cntEdg(thirdType, i) << "; ";
        myfile << "\n";
        std::cout << i << '\n';
    }
}

void start_fordbellman() {
    std::ofstream myfile;
    myfile.open("fordbellman.csv");
    myfile
            << "; full graph; second type graph; third type graph; full graph edg; second type graph edg; third type graph edg;\n";

    for (int i = 10; i <= 1010; i += 50) {
        std::vector<std::vector<int64_t>> fullGraph = generateFullGraph(i);
        std::vector<std::vector<int64_t>> thirdType = generateThirdGraph(i);
        std::vector<std::vector<int64_t>> secondType = generateSecondGraph(i);
        std::vector<std::pair<int64_t, std::pair<int, int>>> fullGraph1;
        for (int i1 = 0; i1 < i; i1++) {
            for (int j = 0; j < i; j++) {
                if (fullGraph[i1][j] != 0) {
                    fullGraph1.push_back({fullGraph[i1][j], {i1, j}});
                }
            }
        }
        std::vector<std::pair<int64_t, std::pair<int, int>>> secondType1;
        for (int i1 = 0; i1 < i; i1++) {
            for (int j = 0; j < i; j++) {
                if (secondType[i1][j] != 0) {
                    secondType1.push_back({secondType[i1][j], {i1, j}});
                }
            }
        }
        std::vector<std::pair<int64_t, std::pair<int, int>>> thirdType1;
        for (int i1 = 0; i1 < i; i1++) {
            for (int j = 0; j < i; j++) {
                if (thirdType[i1][j] != 0) {
                    thirdType1.push_back({thirdType[i1][j], {i1, j}});
                }
            }
        }
        myfile << i << "; ";
        long long ans1 = 0;
        long long ans2 = 0;
        long long ans3 = 0;
        for (int j = 0; j < 5; j++) {
            ans1 += fordbellman(fullGraph1, i);
            ans2 += fordbellman(secondType1, i);
            ans3 += fordbellman(thirdType1, i);
        }
        ans1 /= 5;
        ans2 /= 5;
        ans3 /= 5;
        myfile << ans1 << "; ";
        myfile << ans2 << "; ";
        myfile << ans3 << "; ";
        myfile << cntEdg(fullGraph, i) << "; ";
        myfile << cntEdg(secondType, i) << "; ";
        myfile << cntEdg(thirdType, i) << "; ";
        myfile << "\n";
    }
}

void start_floyd_warshall() {
    std::ofstream myfile;
    myfile.open("floyd_warshall.csv");
    myfile
            << "; full graph; second type graph; third type graph; full graph edg; second type graph edg; third type graph edg; \n";
    //full graph edg; second type graph edg; third type graph edg;
    std::vector<std::vector<int64_t>> fullGraph = generateFullGraph(1010);
    std::vector<std::vector<int64_t>> thirdType = generateThirdGraph(1010);
    for (int i = 10; i <= 1010; i += 50) {
        std::vector<std::vector<int64_t>> secondType = generateSecondGraph(i);
        myfile << i << "; ";
        long long ans1 = 0;
        long long ans2 = 0;
        long long ans3 = 0;
        for (int j = 0; j < 5; j++) {
            ans1 += floyduorshil(fullGraph, i);
            ans2 += floyduorshil(secondType, i);
            ans3 += floyduorshil(thirdType, i);
        }
        ans1 /= 5;
        ans2 /= 5;
        ans3 /= 5;
        myfile << ans1 << "; ";
        myfile << ans2 << "; ";
        myfile << ans3 << "; ";
        myfile << cntEdg(fullGraph, i) << "; ";
        myfile << cntEdg(secondType, i) << "; ";
        myfile << cntEdg(thirdType, i) << "; ";
        myfile << "\n";
    }
}

int main() {
    //start_dejkstra();
    //start_fordbellman();
    //start_floyd_warshall();
    start_johnson();
    return 0;
}
