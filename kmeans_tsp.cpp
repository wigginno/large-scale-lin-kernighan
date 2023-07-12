#include "kmeans_tsp.h"

#include <iostream>
#include <chrono>

// Define a maximum number of iterations for the k-means algorithm
const int MAX_KMEANS_ITERS = 50;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " [problem file]\n";
        return 1;
    }

    // Read problem file
    std::string problem_filename = argv[1];
    std::string solution_filename = problem_filename + ".tour";
    Problem problem(problem_filename, solution_filename);

    auto start = std::chrono::high_resolution_clock::now();
    // Construct the tour using an O(nlogn) insertion heuristic based on k-means clustering
    Tour tour = construct_tour(problem.cities(), MAX_KMEANS_ITERS);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time to construct tour: " << elapsed_ms.count() << " ms\n";

    // Get the tour
    std::vector<id_t> tour_cities = tour.get_tour();

    // Write the tour to the solution file
    problem.write_solution(tour_cities);

    return 0;
}
