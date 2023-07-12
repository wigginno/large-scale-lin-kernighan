#include <cmath>
#include <fstream>
#include <vector>
#include <set>
#include <utility>
#include <random>
#include <map>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <array>
#include <iostream>
#include <stdexcept>
#include <time.h>
#include <sstream>
#include <limits>

const double DOUBLE_MIN = std::numeric_limits<double>::min();

// Type for city indices/identifiers
using id_t = uint32_t;

inline uint64_t interleave(uint32_t x, uint32_t y) {
    // Bit hacking sorcery from
    // https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
    static const uint64_t B[] = {
        0x5555555555555555ULL,
        0x3333333333333333ULL,
        0x0F0F0F0F0F0F0F0FULL,
        0x00FF00FF00FF00FFULL,
        0x0000FFFF0000FFFFULL
    };
    static const unsigned int S[] = {1, 2, 4, 8, 16};

    uint64_t xx = x;
    uint64_t yy = y;

    xx = (xx | (xx << S[4])) & B[4];
    xx = (xx | (xx << S[3])) & B[3];
    xx = (xx | (xx << S[2])) & B[2];
    xx = (xx | (xx << S[1])) & B[1];
    xx = (xx | (xx << S[0])) & B[0];

    yy = (yy | (yy << S[4])) & B[4];
    yy = (yy | (yy << S[3])) & B[3];
    yy = (yy | (yy << S[2])) & B[2];
    yy = (yy | (yy << S[1])) & B[1];
    yy = (yy | (yy << S[0])) & B[0];

    return xx | (yy << 1);
}

uint64_t interleave(float x, float y) {
    // Assume x and y are floats in range [-1., 1.]
    static const int32_t scale_factor = INT32_MAX / 2;

    uint32_t x_int = static_cast<uint32_t>((x + 1.) * scale_factor);
    uint32_t y_int = static_cast<uint32_t>((y + 1.) * scale_factor);

    return interleave(x_int, y_int);
}

uint64_t FNV1a_hash(uint32_t a, uint32_t b, uint32_t c) {
    static const uint64_t FNV_prime = 1099511628211u;
    static const uint64_t FNV_offset_basis = 14695981039346656037u;

    union toBytes {
        uint32_t i;
        unsigned char c[4];
    };

    // Make the hash function commutative (e.g. hash(a, b, c) == hash(b, a, c))
    if (a > b) {
        std::swap(a, b);
    }
    if (b > c) {
        std::swap(b, c);
    }
    if (a > b) {
        std::swap(a, b);
    }

    toBytes input[3];
    input[0].i = a;
    input[1].i = b;
    input[2].i = c;

    uint64_t hash = FNV_offset_basis;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            hash ^= input[i].c[j];
            hash *= FNV_prime;
        }
    }

    return hash;
}

double euclidean_distance_squared(float x1, float y1, float x2, float y2) {
    double x_diff = x1 - x2;
    double y_diff = y1 - y2;
    return x_diff * x_diff + y_diff * y_diff;
}

double euclidean_distance(float x1, float y1, float x2, float y2) {
    return std::sqrt(euclidean_distance_squared(x1, y1, x2, y2));
}

class Point {
 private:
    float _x, _y;

 public:
    Point(float x, float y) : _x(x), _y(y) {}

    float x() const { return _x; }
    float y() const { return _y; }

    float edist(const Point& other) const {
        return euclidean_distance(_x, _y, other._x, other._y);
    }

    double edist_squared(const Point& other) const {
        return edist_squared(other._x, other._y);
    }

    double edist_squared(float x, float y) const {
        return euclidean_distance_squared(_x, _y, x, y);
    }

    bool operator==(const Point& other) const {
        return _x == other._x && _y == other._y;
    }

    bool operator<(const Point& other) const {
        return interleave(_x, _y) < interleave(other._x, other._y);
    }
};

// Type for vector of points
using PointArray = std::vector<Point>;

struct Problem {
 public:
    Problem(const std::string& problem_filename, const std::string& solution_filename)
    : _problem_filename(problem_filename), _solution_filename(solution_filename) {
        std::ifstream file(_problem_filename);
        id_t n;
        file >> n;

        _cities.reserve(n);
        _city_ids.reserve(n);

        double _min_x = std::numeric_limits<double>::max();
        double _max_x = std::numeric_limits<double>::min();
        double _min_y = std::numeric_limits<double>::max();
        double _max_y = std::numeric_limits<double>::min();

        // First pass to validate city ids and get range of coordinates for normalization
        for (id_t i = 0; i < n; ++i) {
            id_t city_id;
            double x, y;
            file >> city_id >> x >> y;

            // Cities should be numbered 0, 1, 2, ..., n - 1
            assert(city_id == i);

            // Update min and max x and y
            _min_x = std::min(_min_x, x);
            _max_x = std::max(_max_x, x);
            _min_y = std::min(_min_y, y);
            _max_y = std::max(_max_y, y);
        }

        file.clear();
        file.seekg(0, std::ios::beg);
        file >> n;

        double x_range = (_max_x - _min_x);
        double y_range = (_max_y - _min_y);
        _scale = std::max(x_range, y_range) * 2.;

        // Second pass to store city ids and normalized coordinates
        for (id_t i = 0; i < n; ++i) {
            id_t city_id;
            double x, y;
            file >> city_id >> x >> y;

            _city_ids.push_back(city_id);

            // Normalize coordinates to range [-1., 1.]
            float x_norm = static_cast<float>((x - _min_x) / _scale - 1.);
            float y_norm = static_cast<float>((y - _min_y) / _scale - 1.);
            _cities.emplace_back(x_norm, y_norm);
        }

        file.close();

        // Sort city ids by Morton code of their coordinates
        std::sort(_city_ids.begin(), _city_ids.end(), [this](id_t a, id_t b) {
            return _cities[a] < _cities[b];
        });

        // City ids were originally 0, 1, 2, ..., n - 1 and are now permuted
        // Sort cities using the same permutation
        std::vector<bool> in_sorted_position(n, false);
        for (id_t i = 0; i < _cities.size(); ++i) {
            if (in_sorted_position[i]) {
                continue;
            }

            id_t prev_j = i;
            id_t j = _city_ids[i];
            while (j != i) {
                std::swap(_cities[prev_j], _cities[j]);
                in_sorted_position[j] = true;
                prev_j = j;
                j = _city_ids[j];
            }

            in_sorted_position[i] = true;
        }
    }

    const PointArray& cities() {
        return _cities;
    }

    double denormalize(float coord) {
        return (coord + 1.) * _scale;
    }

    double compute_tour_weight(const std::vector<id_t>& tour) {
        // Start tour weight with distance from last city to first city
        double x0 = denormalize(_cities[tour[0]].x());
        double y0 = denormalize(_cities[tour[0]].y());
        double x1 = denormalize(_cities[tour[tour.size() - 1]].x());
        double y1 = denormalize(_cities[tour[tour.size() - 1]].y());

        double tour_weight = euclidean_distance(x0, y0, x1, y1);

        // Add distance between each pair of cities in the tour
        for (id_t i = 0; i < tour.size() - 1; ++i) {
            x0 = denormalize(_cities[tour[i]].x());
            y0 = denormalize(_cities[tour[i]].y());
            x1 = denormalize(_cities[tour[i + 1]].x());
            y1 = denormalize(_cities[tour[i + 1]].y());
            tour_weight += euclidean_distance(x0, y0, x1, y1);
        }

        return tour_weight;
    }

    void write_solution(const std::vector<id_t>& tour, double tour_weight = DOUBLE_MIN) {
        // Verify that the tour is valid
        int n_cities = _cities.size();
        assert (tour.size() == n_cities);
        std::vector<bool> visited(n_cities, false);
        for (id_t city : tour) {
            assert(city < n_cities);
            assert(!visited[city]);
            visited[city] = true;
        }

        if (tour_weight == DOUBLE_MIN) {
            tour_weight = compute_tour_weight(tour);
        }

        // Map city ids back to original ids
        std::vector<id_t> tour_corrected_ids;
        tour_corrected_ids.reserve(tour.size());
        for (id_t city_id : tour) {
            tour_corrected_ids.push_back(_city_ids[city_id]);
        }

        std::ofstream output_file(_solution_filename);
        output_file << tour_weight << "\n";

        std::stringstream ss;
        const id_t chunk_size = 100000;
        int counter = 0;

        for (id_t city : tour) {
            ss << city << "\n";
            counter++;

            if (counter == chunk_size) {
                // Write a chunk of cities to the file
                output_file << ss.str();
                ss.str(std::string());  // Clear the stringstream
                ss.clear(); // Clear error flags
                counter = 0;    // Reset counter
            }
        }

        if (counter > 0) {
            // Write remaining cities to the file
            output_file << ss.str();
        }

        output_file.close();
    }

 private:
    const std::string& _problem_filename;
    const std::string& _solution_filename;
    PointArray _cities;
    std::vector<id_t> _city_ids;
    double _scale;
};

class City {
 public:
    City(id_t index) {
        _index = index;
        _neighbor_1 = index;
        _neighbor_2 = index;
        _neighbor_1_weight = 0;
        _neighbor_2_weight = 0;
    }

    id_t index() {
        return _index;
    }

    bool not_in_tour() {
        return _neighbor_1 == _index;
    }

    id_t neighbor_1() {
        return _neighbor_1;
    }

    id_t neighbor_2() {
        return _neighbor_2;
    }

    float neighbor_1_weight() {
        return _neighbor_1_weight;
    }

    float neighbor_2_weight() {
        return _neighbor_2_weight;
    }

    void set_neighbor_1(id_t neighbor_1, float neighbor_1_weight) {
        _neighbor_1 = neighbor_1;
        _neighbor_1_weight = neighbor_1_weight;
    }

    void set_neighbor_2(id_t neighbor_2, float neighbor_2_weight) {
        _neighbor_2 = neighbor_2;
        _neighbor_2_weight = neighbor_2_weight;
    }

    void replace_neighbor(id_t old_neighbor, id_t new_neighbor, float new_neighbor_weight) {
        if (_neighbor_1 == old_neighbor) {
            set_neighbor_1(new_neighbor, new_neighbor_weight);
        } else {
            set_neighbor_2(new_neighbor, new_neighbor_weight);
        }
    }
 private:
    id_t _index;
    id_t _neighbor_1;
    id_t _neighbor_2;
    float _neighbor_1_weight;
    float _neighbor_2_weight;
};

class Tour {
 public:
    Tour(const PointArray& points): _points(points) {
        _cities.reserve(points.size());
        for (id_t i = 0; i < points.size(); ++i) {
            _cities.emplace_back(i);
        }
    }

    void init_subtour(id_t city_1_index, id_t city_2_index) {
        float weight = _points[city_1_index].edist(_points[city_2_index]);
        _cities[city_1_index].set_neighbor_1(city_2_index, weight);
        _cities[city_1_index].set_neighbor_2(city_2_index, weight);
        _cities[city_2_index].set_neighbor_1(city_1_index, weight);
        _cities[city_2_index].set_neighbor_2(city_1_index, weight);
    }

    void insert_adjacent(id_t city_1_index, id_t city_2_index) {
        bool city_1_not_in_tour = _cities[city_1_index].not_in_tour();
        bool city_2_not_in_tour = _cities[city_2_index].not_in_tour();

        if (city_1_not_in_tour == city_2_not_in_tour) {
            return;
        }

        id_t c1 = city_1_not_in_tour ? city_2_index : city_1_index;
        id_t c2 = city_1_not_in_tour ? city_1_index : city_2_index;

        id_t c0 = _cities[c1].neighbor_1();
        float c1_c0_weight = _cities[c1].neighbor_1_weight();
        id_t c3 = _cities[c1].neighbor_2();
        float c1_c3_weight = _cities[c1].neighbor_2_weight();

        float c1_c2_weight = _points[c1].edist(_points[c2]);
        float c2_c0_weight = _points[c2].edist(_points[c0]);
        float c2_c3_weight = _points[c2].edist(_points[c3]);

        float pos1_added_weight = c1_c2_weight + c2_c0_weight - c1_c0_weight;
        float pos2_added_weight = c1_c2_weight + c2_c3_weight - c1_c3_weight;

        if (pos1_added_weight < pos2_added_weight) {
            _cities[c1].set_neighbor_1(c2, c1_c2_weight);
            _cities[c0].replace_neighbor(c1, c2, c2_c0_weight);
            _cities[c2].set_neighbor_1(c0, c2_c0_weight);
            _cities[c2].set_neighbor_2(c1, c1_c2_weight);
        } else {
            _cities[c1].set_neighbor_2(c2, c1_c2_weight);
            _cities[c3].replace_neighbor(c1, c2, c2_c3_weight);
            _cities[c2].set_neighbor_1(c1, c1_c2_weight);
            _cities[c2].set_neighbor_2(c3, c2_c3_weight);
        }
    }

    std::vector<id_t> get_tour() {
        std::vector<id_t> tour;
        tour.reserve(_cities.size());
        std::vector<bool> visited(_cities.size(), false);
        tour.push_back(0);
        visited[0] = true;
        id_t next_point_index = _get_next_point_index(0, visited);

        while (!visited[next_point_index]) {
            tour.push_back(next_point_index);
            visited[next_point_index] = true;
            next_point_index = _get_next_point_index(next_point_index, visited);
        }

        return tour;
    }

 private:
    std::vector<City> _cities;
    const PointArray& _points;

    id_t _get_next_point_index(id_t point_index, std::vector<bool>& visited) {
        id_t next_point_index = _cities[point_index].neighbor_1();
        if (visited[next_point_index]) {
            next_point_index = _cities[point_index].neighbor_2();
        }
        return next_point_index;
    }
};

class Cluster {
 private:
    std::vector<id_t> _point_indices;
    id_t _centroid;

 public:
    Cluster(id_t centroid) : _centroid(centroid) {}

    id_t centroid() {
        return _centroid;
    }

    void update_centroid(const PointArray& points) {
        double x = 0.0;
        double y = 0.0;

        for (id_t i = 0; i < size(); ++i) {
            const Point& point = points[_point_indices[i]];
            x += point.x();
            y += point.y();
        }

        x /= size();
        y /= size();

        // return point with min distance to center of cluster
        double min_dist = std::numeric_limits<double>::max();
        id_t min_index = 0;
        for (id_t i = 0; i < size(); ++i) {
            const Point& point = points[_point_indices[i]];
            double dist = point.edist_squared(x, y);
            if (dist < min_dist) {
                min_dist = dist;
                min_index = i;
            }
        }

        _centroid = _point_indices[min_index];
    }

    size_t size() const {
        return _point_indices.size();
    }

    std::vector<uint32_t>::iterator begin() {
        return _point_indices.begin();
    }

    std::vector<uint32_t>::iterator end() {
        return _point_indices.end();
    }

    void add_point(uint32_t point) {
        _point_indices.push_back(point);
    }

    void clear_points() {
        _point_indices.clear();
    }

    uint32_t rand_point() {
        size_t index = std::rand() % size();
        return _point_indices[index];
    }

    uint32_t farthest_point(uint32_t point_index, const PointArray& points) {
        double max_dist = DOUBLE_MIN;
        uint32_t max_index = 0;
        const Point& point = points[point_index];
        for (id_t i = 0; i < size(); ++i) {
            double dist = point.edist_squared(points[_point_indices[i]]);
            if (dist > max_dist) {
                max_dist = dist;
                max_index = i;
            }
        }

        return _point_indices[max_index];
    }
};

std::pair<Cluster, Cluster> divide_cluster(Cluster& original_cluster, const PointArray& points, int max_kmeans_iters) {
    uint32_t rand_point_index = original_cluster.rand_point();

    Cluster cluster_1(rand_point_index);
    Cluster cluster_2(original_cluster.farthest_point(rand_point_index, points));

    for (int i = 0; i < max_kmeans_iters; ++i) {
        id_t centroid_1_index = cluster_1.centroid();
        id_t centroid_2_index = cluster_2.centroid();
        const Point& centroid_1 = points[centroid_1_index];
        const Point& centroid_2 = points[centroid_2_index];

        for (const auto& point_index : original_cluster) {
            const Point& point = points[point_index];
            if (point_index == centroid_1_index) {
                cluster_1.add_point(point_index);
            } else if (point_index == centroid_2_index) {
                cluster_2.add_point(point_index);
            } else if (point.edist_squared(centroid_1) < point.edist_squared(centroid_2)) {
                cluster_1.add_point(point_index);
            } else {
                cluster_2.add_point(point_index);
            }
        }

        cluster_1.update_centroid(points);
        cluster_2.update_centroid(points);

        if (cluster_1.centroid() == centroid_1_index && cluster_2.centroid() == centroid_2_index) {
            // No change in centroids, so we are done
            break;
        }

        if (i != max_kmeans_iters - 1) {
            cluster_1.clear_points();
            cluster_2.clear_points();
        }
    }

    return std::make_pair(std::move(cluster_1), std::move(cluster_2));
}

void construct(Tour& tour, Cluster& current_cluster, const PointArray& points, int max_kmeans_iters) {
    auto subclusters = divide_cluster(current_cluster, points, max_kmeans_iters);
    id_t centroid_0 = current_cluster.centroid();
    id_t centroid_1 = subclusters.first.centroid();
    id_t centroid_2 = subclusters.second.centroid();
    tour.insert_adjacent(centroid_0, centroid_1);
    tour.insert_adjacent(centroid_0, centroid_2);

    if (subclusters.first.size() > 1) {
        construct(tour, subclusters.first, points, max_kmeans_iters);
    }

    if (subclusters.second.size() > 1) {
        construct(tour, subclusters.second, points, max_kmeans_iters);
    }
}

Tour construct_tour(const PointArray& points, int max_kmeans_iters) {
    srand(time(NULL));
    Tour tour(points);
    Cluster cluster(0);

    for (id_t i = 0; i < points.size(); ++i) {
        cluster.add_point(i);
    }

    auto clusters = divide_cluster(cluster, points, max_kmeans_iters);
    Cluster& cluster_1 = clusters.first;
    Cluster& cluster_2 = clusters.second;

    tour.init_subtour(cluster_1.centroid(), cluster_2.centroid());

    if (cluster_1.size() > 1) {
        construct(tour, cluster_1, points, max_kmeans_iters);
    }

    if (cluster_2.size() > 1) {
        construct(tour, cluster_2, points, max_kmeans_iters);
    }

    return tour;
}
