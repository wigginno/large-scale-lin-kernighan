#include "delaunay.h"

#include <iostream>
#include <chrono>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " [problem file]\n";
        return 1;
    }

    // Read problem file
    std::string problem_filename = argv[1];
    std::string solution_filename = problem_filename + ".tri";
    auto t1 = std::chrono::high_resolution_clock::now();
    Problem problem(problem_filename, solution_filename);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "Time to read problem file: " << t_elapsed.count() << " ms\n";

    auto start = std::chrono::high_resolution_clock::now();
    // Construct the Delaunay triangulation using randomized incremental insertion
    std::vector<std::tuple<id_t, id_t, id_t>> simplices;
    Triangulation triangulation(problem.cities(), simplices);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Total time to construct Delaunay triangulation: " << elapsed_ms.count() << " ms\n";

    // Write the triangulation to the solution file
    problem.write_delaunay(simplices);
    return 0;
}
