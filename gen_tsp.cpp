#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sstream>

int main(int argc, char *argv[]) {
    if(argc != 3) {
        std::cerr << "Usage: gen_tsp n output_file" << std::endl;
        return 1;
    }

    size_t n = std::stoi(argv[1]);
    std::string output_file = argv[2];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, n);

    std::set<std::pair<uint32_t, uint32_t>> points;

    while (points.size() < n) {
        uint32_t x = distrib(gen);
        uint32_t y = distrib(gen);
        points.insert({x, y});
    }

    std::ofstream out(output_file);
    if (!out) {
        std::cerr << "Could not open output file." << std::endl;
        return 1;
    }

    // get output filename without prefix or extension
    std::string filename = output_file.substr(output_file.find_last_of("/\\") + 1);
    filename = filename.substr(0, filename.find_last_of("."));
    out << "NAME: " << filename << '\n';
    out << "TYPE: TSP\n";
    out << "DIMENSION: " << n << '\n';
    out << "EDGE_WEIGHT_TYPE: EUC_2D\n";
    out << "NODE_COORD_SECTION\n";
    int idx = 1;
    for (const auto& point : points) {
        out << idx++ << ' ' << point.first << ' ' << point.second << '\n';
    }

    return 0;
}
