#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <chrono>

// Type for point indices/identifiers.
using id_t = uint32_t;

// Type for vector of point indices/identifiers.
using id_vec = std::vector<id_t>;

inline uint64_t interleave(uint32_t x, uint32_t y) {
    // Wizardry from https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
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

class Point {
 private:
    float _x, _y;

 public:
    Point(float x, float y) : _x(x), _y(y) {}

    float x() const { return _x; }
    float y() const { return _y; }

    bool operator<(const Point& other) const {
        return interleave(_x, _y) < interleave(other._x, other._y);
    }
};

// Type for vector of points
using PointArray = std::vector<Point>;

/**
 * Hash function for unordered_map with keys that are already uint64 hashes or unique ids.
*/
template<typename K>
struct NoHash {
    K operator()(const K& key) const {
        return key;
    }
};

// Type for unordered_set of unique ids.
using id_unordered_set = std::unordered_set<id_t, NoHash<id_t>>;

class Triangle {
 public:
    Triangle(std::tuple<id_t, id_t, id_t> point_ids)
        : a_id(std::get<0>(point_ids)), b_id(std::get<1>(point_ids)), c_id(std::get<2>(point_ids)) {}

    Triangle(id_t a_id, id_t b_id, id_t c_id)
        : a_id(a_id), b_id(b_id), c_id(c_id) {}

    const id_t a_id, b_id, c_id;

    bool contains(id_t point_id) const {
        return point_id == a_id || point_id == b_id || point_id == c_id;
    }

    id_vec intersection(const Triangle& other) const {
        return other.contains(a_id) 
            ? (other.contains(b_id) 
                ? (other.contains(c_id) ? id_vec{a_id, b_id, c_id} : id_vec{a_id, b_id})
                : (other.contains(c_id) ? id_vec{a_id, c_id} : id_vec{a_id}))
            : (other.contains(b_id)
                ? (other.contains(c_id) ? id_vec{b_id, c_id} : id_vec{b_id})
                : (other.contains(c_id) ? id_vec{c_id} : id_vec{}));
    }
};

uint64_t FNV1a_hash(id_t a, id_t b, id_t c) {
    static const uint64_t FNV_prime = 1099511628211u;
    static const uint64_t FNV_offset_basis = 14695981039346656037u;

    union toBytes {
        id_t i;
        unsigned char c[4];
    };

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

uint64_t commutative_FNV1a_hash(id_t a, id_t b, id_t c) {
    if (a > c) std::swap(a, c);
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    return FNV1a_hash(a, b, c);
}

uint64_t triangle_hash(const Triangle& triangle) {
    return commutative_FNV1a_hash(triangle.a_id, triangle.b_id, triangle.c_id);
}

uint64_t triangle_hash(const Triangle* triangle) {
    return commutative_FNV1a_hash(triangle->a_id, triangle->b_id, triangle->c_id);
}

uint64_t triangle_hash(id_t a_id, id_t b_id, id_t c_id) {
    return commutative_FNV1a_hash(a_id, b_id, c_id);
}

enum Orientation {
    CCW = 1,
    CW = 2,
    COLLINEAR = 3
};

class TriPoint {
 public:
    id_t id;
    float x, y;

    TriPoint(id_t id, float x, float y)
        : id(id), x(x), y(y) {}

    TriPoint(id_t id, const Point& point)
        : id(id), x(point.x()), y(point.y()) {}

    void add_neighbor(id_t neighbor_id) {
        _neighbors.push_back(neighbor_id);
    }

    void remove_neighbor(id_t neighbor_id) {
        _neighbors.erase(
            std::remove(_neighbors.begin(), _neighbors.end(), neighbor_id), _neighbors.end());
    }

    void add_neighbors(id_t neighbor_id_1, id_t neighbor_id_2) {
        add_neighbor(neighbor_id_1);
        add_neighbor(neighbor_id_2);
    }

    void add_neighbors(id_t neighbor_id_1, id_t neighbor_id_2, id_t neighbor_id_3) {
        add_neighbor(neighbor_id_1);
        add_neighbor(neighbor_id_2);
        add_neighbor(neighbor_id_3);
    }

    bool has_neighbor(id_t neighbor_id) const {
        return std::find(_neighbors.begin(), _neighbors.end(), neighbor_id) != _neighbors.end();
    }

    id_vec intersection(const TriPoint* other) {
        id_vec intersection;
        intersection.reserve(_neighbors.size());
        for (id_t neighbor_id : _neighbors) {
            if (other->has_neighbor(neighbor_id)) {
                intersection.push_back(neighbor_id);
            }
        }
        return intersection;
    }

 private:
    id_vec _neighbors;
};

Orientation orientation(const TriPoint* A, const TriPoint* B, const TriPoint* C) {
    double det = (static_cast<double>(A->x) - C->x) * (static_cast<double>(B->y) - C->y)
                 - (static_cast<double>(A->y) - C->y) * (static_cast<double>(B->x) - C->x);

    if (std::abs(det) < 1e-9) {
        return COLLINEAR;
    }

    return det > 0 ? CCW : CW;
}

struct InTriangleTest {
    const double A_x, A_y, B_x, B_y, C_x, C_y, T_1, T_2, T_3;
    InTriangleTest(const TriPoint* A, const TriPoint* B, const TriPoint* C)
        : A_x(static_cast<double>(A->x)), A_y(static_cast<double>(A->y)),
          B_x(static_cast<double>(B->x)), B_y(static_cast<double>(B->y)),
          C_x(static_cast<double>(C->x)), C_y(static_cast<double>(C->y)),
          T_1(A_y * B_x - A_x * B_y - 1e-9),
          T_2(B_y * C_x - B_x * C_y - 1e-9),
          T_3(C_y * A_x - C_x * A_y - 1e-9) {}

    bool operator()(float D_x, float D_y) const {
        double T_4 = A_x * D_y - A_y * D_x;
        double T_5 = B_x * D_y - B_y * D_x;
        double T_6 = C_x * D_y - C_y * D_x;

        return (T_5 - T_4 > T_1 && T_6 - T_5 > T_2 && T_4 - T_6 > T_3);
    }
};

bool is_in_circle(const TriPoint* A, const TriPoint* B, const TriPoint* C, const TriPoint* D) {
    double AD_x = static_cast<double>(A->x) - D->x;
    double AD_y = static_cast<double>(A->y) - D->y;
    double BD_x = static_cast<double>(B->x) - D->x;
    double BD_y = static_cast<double>(B->y) - D->y;
    double CD_x = static_cast<double>(C->x) - D->x;
    double CD_y = static_cast<double>(C->y) - D->y;

    double det = (
        (AD_x * AD_x + AD_y * AD_y) * (BD_x * CD_y - BD_y * CD_x) -
        (BD_x * BD_x + BD_y * BD_y) * (AD_x * CD_y - AD_y * CD_x) +
        (CD_x * CD_x + CD_y * CD_y) * (AD_x * BD_y - AD_y * BD_x)
    );

    return det > 1e-9;
}

class Triangulation {
 public:
    Triangulation(const PointArray& point_coords, std::vector<std::tuple<id_t, id_t, id_t>>& simplices)
        : _n_vertices(point_coords.size()), _super_triangle(_get_super_triangle_ids())
    {
        _init_triangulation(point_coords);
        _init_shuffled_point_ids();
        _construct();
        _get_simplices(simplices);
    }
 private:
    const id_t _n_vertices;
    Triangle _super_triangle;
    std::vector<id_t> _shuffled_point_ids;
    std::vector<std::unique_ptr<TriPoint>> _vertices;
    std::unordered_map<uint64_t, Triangle, NoHash<uint64_t>> _triangles;
    std::unordered_map<uint64_t, std::unique_ptr<id_unordered_set>, NoHash<uint64_t>> _buckets;     // TODO: should point to vectors not sets
    std::vector<const Triangle*> _containing_triangle;

    // TODO: Make all current set/map data members vectors except for one unordered map of triangle hash (uint64 key) to index (id_t value)

    void _init_shuffled_point_ids() {
        _shuffled_point_ids.resize(_n_vertices);
        std::iota(_shuffled_point_ids.begin(), _shuffled_point_ids.end(), 0);
        std::mt19937 g(12345);
        std::shuffle(_shuffled_point_ids.begin(), _shuffled_point_ids.end(), g);
    }

    void _init_triangulation(const PointArray& point_coords) {
        std::unique_ptr<id_unordered_set> initial_bucket = std::make_unique<id_unordered_set>();
        initial_bucket->reserve(_n_vertices);
        for (id_t i = 0; i < _n_vertices; i++) {
            initial_bucket->insert(i);
        }

        _buckets[triangle_hash(_super_triangle.a_id, _super_triangle.b_id, _super_triangle.c_id)]
            = std::move(initial_bucket);

        _vertices.reserve(_n_vertices + 3);
        _containing_triangle.reserve(_n_vertices);

        for (id_t i = 0; i < _n_vertices; i++) {
            _vertices.emplace_back(std::make_unique<TriPoint>(i, point_coords[i]));
            _containing_triangle.push_back(&_super_triangle);
        }

        id_t a_id = _super_triangle.a_id;
        id_t b_id = _super_triangle.b_id;
        id_t c_id = _super_triangle.c_id;
        _vertices.emplace_back(std::make_unique<TriPoint>(a_id, -2.43, -1.01));
        _vertices.emplace_back(std::make_unique<TriPoint>(b_id, 2.43, -1.01));
        _vertices.emplace_back(std::make_unique<TriPoint>(c_id, 0, 2.43));
        _vertices[a_id]->add_neighbors(b_id, c_id);
        _vertices[b_id]->add_neighbors(a_id, c_id);
        _vertices[c_id]->add_neighbors(a_id, b_id);

        _triangles.emplace(triangle_hash(a_id, b_id, c_id), _super_triangle);
    }

    const TriPoint* _get_point(id_t id) const {
        return _vertices[id].get();
    }

    std::tuple<id_t, id_t, id_t> _get_super_triangle_ids() const {
        return std::make_tuple(_n_vertices, _n_vertices + 1, _n_vertices + 2);
    }

    std::tuple<id_t, id_t, id_t> _get_ccw_triangle_ids(id_t a_id, id_t b_id, id_t c_id) const {
        return orientation(_get_point(a_id), _get_point(b_id), _get_point(c_id)) == CCW
            ? std::make_tuple(a_id, b_id, c_id)
            : (orientation(_get_point(b_id), _get_point(a_id), _get_point(c_id)) == CCW
                ? std::make_tuple(b_id, a_id, c_id)
                : std::make_tuple(c_id, a_id, b_id));
    }

    id_t _get_opposite_vertex_id(id_t a_id, id_t b_id, id_t p_id) const {
        id_vec intersection = _vertices[a_id]->intersection(_get_point(b_id));
        for (id_t id : intersection) {
            if (id != p_id && _triangles.count(triangle_hash(a_id, b_id, id)) > 0) {
                return id;
            }
        }
        throw std::runtime_error("Opposite vertex not found");
    }

    void _update_buckets(const std::vector<Triangle*>& new_triangles,
                         const std::vector<uint64_t>& old_triangle_hashes,
                         const std::vector<uint64_t>& new_triangle_hashes) {
        // Add new buckets
        for (uint64_t new_triangle_hash : new_triangle_hashes) {
            _buckets.emplace(new_triangle_hash, std::make_unique<id_unordered_set>());
        }

        id_t last_new_triangle_index = new_triangles.size() - 1;

        // Re-bucket points from old triangles to new triangles
        for (uint64_t old_triangle_hash : old_triangle_hashes) {
            //std::vector<const TriPoint*> points;
            //points.reserve(_buckets[old_triangle_hash]->size());
            //for (id_t point_id : *_buckets[old_triangle_hash]) {
            //    points.push_back(_get_point(point_id));
            //}
            std::vector<const TriPoint*> points(_buckets[old_triangle_hash]->size());
            std::transform(std::execution::par, _buckets[old_triangle_hash]->begin(), _buckets[old_triangle_hash]->end(),
                points.begin(),
                [&](id_t point_id) {
                    return _get_point(point_id);
                }
            );

            _buckets.erase(old_triangle_hash);

            std::vector<std::unique_ptr<InTriangleTest>> tests(new_triangles.size());

            std::transform(std::execution::par, new_triangles.begin(), new_triangles.begin() + last_new_triangle_index,
                tests.begin(),
                [&](Triangle* tri) {
                    return std::make_unique<InTriangleTest>(_get_point(tri->a_id), _get_point(tri->b_id), _get_point(tri->c_id));
                }
            );

            std::for_each(std::execution::par, points.begin(), points.end(),
                [&](const TriPoint* point) {
                    for (id_t i = 0; i < new_triangles.size(); i++) {
                        if (i == last_new_triangle_index || tests[i]->operator()(point->x, point->y)) {
                            _buckets[new_triangle_hashes[i]]->insert(point->id);
                            _containing_triangle[point->id] = new_triangles[i];
                            break;
                        }
                    }
                }
            );
        }
    }

    /**
     * Replace edge ab with edge pd.
    **/
    void _flip_edge(id_t a_id, id_t b_id, id_t p_id, id_t d_id) {
        // Get triangle hashes
        std::vector<uint64_t> old_triangle_hashes = {
            triangle_hash(a_id, b_id, p_id),
            triangle_hash(a_id, b_id, d_id)
        };
        std::vector<uint64_t> new_triangle_hashes = {
            triangle_hash(p_id, d_id, a_id),
            triangle_hash(p_id, d_id, b_id)
        };

        // Remove old triangles
        _triangles.erase(old_triangle_hashes[0]);
        _triangles.erase(old_triangle_hashes[1]);

        // Update neighbors
        _vertices[a_id]->remove_neighbor(b_id);
        _vertices[b_id]->remove_neighbor(a_id);
        _vertices[p_id]->add_neighbor(d_id);
        _vertices[d_id]->add_neighbor(p_id);

        // Add new triangles and prepare vector of pointers to them
        std::vector<Triangle*> new_triangles = {
            &_triangles.emplace(new_triangle_hashes[0], Triangle(_get_ccw_triangle_ids(p_id, d_id, a_id))).first->second,
            &_triangles.emplace(new_triangle_hashes[1], Triangle(_get_ccw_triangle_ids(p_id, d_id, b_id))).first->second
        };

        // Update buckets
        _update_buckets(new_triangles, old_triangle_hashes, new_triangle_hashes);
    }

    void _swap_test(id_t a_id, id_t b_id, id_t p_id) {
        if (a_id >= _n_vertices && b_id >= _n_vertices) {
            return;
        }
        id_t d_id = _get_opposite_vertex_id(a_id, b_id, p_id);
        if (is_in_circle(_get_point(b_id), _get_point(p_id), _get_point(a_id), _get_point(d_id))) {
            _flip_edge(a_id, b_id, p_id, d_id);
            _swap_test(a_id, d_id, p_id);
            _swap_test(d_id, b_id, p_id);
        }
    }

    void _insert_point(id_t p_id) {
        // Find containing triangle
        const Triangle* abc = _containing_triangle[p_id];
        id_t a_id = abc->a_id;
        id_t b_id = abc->b_id;
        id_t c_id = abc->c_id;
        uint64_t abc_hash = triangle_hash(a_id, b_id, c_id);

        _buckets[abc_hash]->erase(p_id);
        _triangles.erase(abc_hash);

        // Update neighbors
        _vertices[a_id]->add_neighbor(p_id);
        _vertices[b_id]->add_neighbor(p_id);
        _vertices[c_id]->add_neighbor(p_id);
        _vertices[p_id]->add_neighbors(a_id, b_id, c_id);

        // Get triangle hashes
        std::vector<uint64_t> old_triangle_hashes = {abc_hash};
        std::vector<uint64_t> new_triangle_hashes = {
            triangle_hash(p_id, a_id, b_id),
            triangle_hash(p_id, a_id, c_id),
            triangle_hash(p_id, b_id, c_id)
        };

        // Add new triangles and prepare vector of pointers to them
        std::vector<Triangle*> new_triangles = {
            &_triangles.emplace(new_triangle_hashes[0], Triangle(_get_ccw_triangle_ids(p_id, a_id, b_id))).first->second,
            &_triangles.emplace(new_triangle_hashes[1], Triangle(_get_ccw_triangle_ids(p_id, a_id, c_id))).first->second,
            &_triangles.emplace(new_triangle_hashes[2], Triangle(_get_ccw_triangle_ids(p_id, b_id, c_id))).first->second
        };

        // Update buckets
        _update_buckets(new_triangles, old_triangle_hashes, new_triangle_hashes);

        // Fix triangulation
        _swap_test(a_id, b_id, p_id);
        _swap_test(b_id, c_id, p_id);
        _swap_test(c_id, a_id, p_id);
    }

    void _construct() {
        for (id_t p_id : _shuffled_point_ids) {
            // Insert point
            _insert_point(p_id);
        }

        // Remove super triangle
        for (auto it = _triangles.begin(); it != _triangles.end();) {
            if (it->second.intersection(_super_triangle).size() > 0) {
                it = _triangles.erase(it);
            } else {
                ++it;
            }
        }
    }

    void _get_simplices(std::vector<std::tuple<id_t, id_t, id_t>>& simplices) const {
        simplices.reserve(_triangles.size());
        for (const auto& triangle_kv : _triangles) {
            const Triangle& triangle = triangle_kv.second;
            simplices.emplace_back(triangle.a_id, triangle.b_id, triangle.c_id);
        }
    }
};

struct Problem {
 public:
    Problem(const std::string& problem_filename, const std::string& solution_filename)
    : _problem_filename(problem_filename), _solution_filename(solution_filename) {
        std::ifstream file(_problem_filename);

        _read_first_city_id(file);

        id_t n = _read_header(file);
        _cities.reserve(n);
        _city_ids.reserve(n);

        double _min_x = std::numeric_limits<double>::max();
        double _max_x = std::numeric_limits<double>::min();
        double _min_y = std::numeric_limits<double>::max();
        double _max_y = std::numeric_limits<double>::min();

        std::string line;
        // First pass to get range of coordinates for normalization
        while (std::getline(file, line)) {
            std::istringstream is(line);
            id_t city_id;
            double x, y;

            if (is >> city_id >> x >> y) {
                // Update min and max x and y
                _min_x = std::min(_min_x, x);
                _max_x = std::max(_max_x, x);
                _min_y = std::min(_min_y, y);
                _max_y = std::max(_max_y, y);
            }
        }

        file.clear();
        file.seekg(0, std::ios::beg);
        _skip_header(file);

        double x_range = _max_x - _min_x;
        double y_range = _max_y - _min_y;
        _scale = std::max(x_range, y_range);

        // Second pass to store city ids and normalized coordinates
        while (std::getline(file, line)) {
            std::istringstream is(line);
            id_t city_id;
            double x, y;

            if (is >> city_id >> x >> y) {
                _city_ids.push_back(city_id - _first_city_id);
                // Normalize coordinates to range [-1., 1.]
                float x_norm = static_cast<float>(((x - _min_x) / _scale) * 2 - 1);
                float y_norm = static_cast<float>(((y - _min_y) / _scale) * 2 - 1);
                _cities.emplace_back(x_norm, y_norm);
            }
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

    void write_delaunay(const std::vector<std::tuple<id_t, id_t, id_t>>& simplices) {
        std::vector<std::tuple<id_t, id_t, id_t>> simplices_mapped_ids;
        simplices_mapped_ids.reserve(simplices.size());
        for (const auto& simplex : simplices) {
            simplices_mapped_ids.emplace_back(
                _city_ids[std::get<0>(simplex)],
                _city_ids[std::get<1>(simplex)],
                _city_ids[std::get<2>(simplex)]
            );
        }

        std::ofstream output_file(_solution_filename);
        std::stringstream ss;
        const id_t chunk_size = 100000;
        int counter = 0;

        for (const auto& simplex : simplices_mapped_ids) {
            ss << std::get<0>(simplex) << " " << std::get<1>(simplex) << " " << std::get<2>(simplex) << "\n";
            counter++;
            if (counter == chunk_size) {
                // Write a chunk of simplices to the file
                output_file << ss.str();
                ss.str(std::string());  // Clear the stringstream
                ss.clear(); // Clear error flags
                counter = 0;    // Reset counter
            }
        }

        if (counter > 0) {
            // Write remaining simplices to the file
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
    id_t _first_city_id;

    id_t _read_header(std::ifstream& file) {
        std::string line;
        id_t n = 0;
        while (std::getline(file, line)) {
            std::string temp;
            std::istringstream is(line);
            if (is >> temp && (temp == "DIMENSION" || temp == "dimension"
                               || temp == "DIMENSION:" || temp == "dimension:")) {
                if (is.peek() == ':') {
                    is.ignore();  // ignore the colon
                }
                if (is >> n) {
                    // Successfully read the number of cities
                } else {
                    throw std::runtime_error("Invalid DIMENSION line in file");
                }
            }
            if (line == "NODE_COORD_SECTION") {
                break;
            }
        }
        return n;
    }

    void _skip_header(std::ifstream& file) {
        std::string line;
        while (std::getline(file, line)) {
            if (line == "NODE_COORD_SECTION") {
                break;
            }
        }
    }

    void _reset_file(std::ifstream& file) {
        file.clear();
        file.seekg(0, std::ios::beg);
    }

    // Read first city id from file and go back to the beginning of the file.
    void _read_first_city_id(std::ifstream& file) {
        _skip_header(file);
        std::string line;
        std::getline(file, line);
        std::istringstream is(line);
        is >> _first_city_id;
        _reset_file(file);
    }
};
