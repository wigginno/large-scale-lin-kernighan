/*
disclaimer: this is an unreadable mess written as a prototype

pros:
- faster than other serial insertion algorithms (~5 seconds for 1M points)
- easy to adapt for solving massive point sets using streaming/out-of-core/partitioning techniques
- potential for further algorithmic improvements
cons:
- tricky to parallelize effectively
- ~5x slower than divide-and-conquer algorithms

compile:
g++ -O3 -march=native -std=c++17 -o delaunay delaunay.cpp
debug:
g++ -g -std=c++17 -o delaunay delaunay.cpp
*/

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <random>

double jitter_if_zero(double value);

struct PerpendicularBisector;

struct Point {
    double x;
    double y;
    Point midpoint(const Point& other_point) const;
    double distance_squared(const Point& other_point) const;
    double distance(const Point& other_point) const;
    double orientation(const Point& a, const Point& b) const;
    Point();
    Point(double x, double y);
    Point(const Point& b, const Point& c, const PerpendicularBisector& ab_perp);
    double operator[](uint32_t i) const {
        return i == 0 ? x : y;
    }
};

struct Line {
    double slope;
    double intercept;

    Line() : slope(0), intercept(0) {}
    Line(const Point& a, const Point& b) {
        slope = jitter_if_zero(b.y - a.y) / jitter_if_zero(b.x - a.x);
        intercept = a.y - slope * a.x;
    }
};

struct PerpendicularBisector : Line {
    void compute_perpendicular_bisector(double diff_x, double diff_y, const Point& midpoint);
    PerpendicularBisector();
    PerpendicularBisector(const Point& a, const Point& b);
};

struct BBox {
    double x_min, x_max, y_min, y_max;

    BBox(double x_min, double x_max, double y_min, double y_max)
    : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max) {}

    BBox() : x_min(-1), x_max(1), y_min(-1), y_max(1) {}

    // Constructor for updating boundaries
    BBox(const BBox& other, double val, int axis, bool is_right) {
        if (is_right) {
            set_right(other, val, axis);
        } else {
            set_left(other, val, axis);
        }
    }

    void set_left(const BBox& other, double val, int axis) {
        if (axis == 0) {
            x_min = other.x_min;
            x_max = val;
            y_min = other.y_min;
            y_max = other.y_max;
        } else {
            x_min = other.x_min;
            x_max = other.x_max;
            y_min = other.y_min;
            y_max = val;
        }
    }

    void set_right(const BBox& other, double val, int axis) {
        if (axis == 0) {
            x_min = val;
            x_max = other.x_max;
            y_min = other.y_min;
            y_max = other.y_max;
        } else {
            x_min = other.x_min;
            x_max = other.x_max;
            y_min = val;
            y_max = other.y_max;
        }
    }

    void reset() {
        x_min = -1;
        x_max = 1;
        y_min = -1;
        y_max = 1;
    }

    bool x_in_bound(double x) const {
        return x >= x_min && x <= x_max;
    }

    bool y_in_bound(double y) const {
        return y >= y_min && y <= y_max;
    }

    // Method for line and bounding box intersection check
    bool intersects(const Line& line) const {
        return y_in_bound(line.slope * x_min + line.intercept)
            || y_in_bound(line.slope * x_max + line.intercept)
            || x_in_bound((y_min - line.intercept) / line.slope)
            || x_in_bound((y_max - line.intercept) / line.slope);
    }

    double dist_to_point(const Point& p) const {
        if (p.x < x_min) {
            if (p.y < y_min) {
                return pow(p.x - x_min, 2) + pow(p.y - y_min, 2);
            } else if (p.y > y_max) {
                return pow(p.x - x_min, 2) + pow(p.y - y_max, 2);
            }
            return pow(p.x - x_min, 2);
        }
        if (p.x > x_max) {
            if (p.y < y_min) {
                return pow(p.x - x_max, 2) + pow(p.y - y_min, 2);
            } else if (p.y > y_max) {
                return pow(p.x - x_max, 2) + pow(p.y - y_max, 2);
            }
            return pow(p.x - x_max, 2);
        }
        if (p.y < y_min) {
            return pow(p.y - y_min, 2);
        } else if (p.y > y_max) {
            return pow(p.y - y_max, 2);
        }
        return 0;
    }
};

void PerpendicularBisector::compute_perpendicular_bisector(double diff_x, double diff_y, const Point& midpoint) {
    slope = diff_x / diff_y;
    intercept = midpoint.y - slope * midpoint.x;
}

struct PerpendicularUnitVector {
    double x;
    double y;
    PerpendicularBisector line;
    PerpendicularUnitVector(const Point& a, const Point& b, int orientation_sign);
};

struct KDNode {
    uint32_t left_index;
    uint32_t right_index;
    void set_left(uint32_t idx) { left_index = idx; }
    void set_right(uint32_t idx) { right_index = idx; }
};

struct KDTree {
    std::vector<KDNode> nodes;
    uint32_t root;
    uint32_t max_depth;
    KDTree(const std::vector<Point>& points);
};

Point::Point() {}

Point::Point(double x, double y) : x(x), y(y) {}

Point Point::midpoint(const Point& other_point) const {
    return {(x + other_point.x) / 2, (y + other_point.y) / 2};
}

double Point::orientation(const Point& a, const Point& b) const {
    return (a.x - x) * (b.y - y) - (a.y - y) * (b.x - x);
}

double Point::distance_squared(const Point& other_point) const {
    return pow(x - other_point.x, 2) + pow(y - other_point.y, 2);
}

double Point::distance(const Point& other_point) const {
    return sqrt(distance_squared(other_point));
}

// Return random jitter between 1e-12 and 5e-12 if (abs(value) < 10^-12).
double jitter_if_zero(double value) {
    while (abs(value) < 1e-12) {
        value = (static_cast<double>(rand()) / RAND_MAX - 0.5) * 1e-11;
    }
    return value;
}

PerpendicularBisector::PerpendicularBisector() {}

PerpendicularBisector::PerpendicularBisector(const Point& a, const Point& b) {
    double diff_x = jitter_if_zero(b.x - a.x);
    double diff_y = jitter_if_zero(a.y - b.y);
    compute_perpendicular_bisector(diff_x, diff_y, a.midpoint(b));
}

// Special `Point` constructor to compute the center of a circumscribed circle of a triangle.
Point::Point(const Point& b, const Point& c, const PerpendicularBisector& ab_perp) {
    PerpendicularBisector bc_perp(b, c);
    double bisector_slope_diff = jitter_if_zero(ab_perp.slope - bc_perp.slope);
    x = (bc_perp.intercept - ab_perp.intercept) / bisector_slope_diff;
    y = x * ab_perp.slope + ab_perp.intercept;
}

PerpendicularUnitVector::PerpendicularUnitVector(const Point& a, const Point& b, int orientation_sign) {
    double diff_x = jitter_if_zero(b.x - a.x);
    double diff_y = jitter_if_zero(a.y - b.y);
    line.compute_perpendicular_bisector(diff_x, diff_y, a.midpoint(b));
    double distance = sqrt(diff_x * diff_x + diff_y * diff_y);
    x = orientation_sign * diff_y / distance;
    y = orientation_sign * diff_x / distance;
}

bool compare_points_for_graham_scan(const Point& a, const Point& b, const Point& p) {
    double ap_dx = a.x - p.x;
    double ap_dy = a.y - p.y;
    double bp_dx = b.x - p.x;
    double bp_dy = b.y - p.y;

    if (ap_dy == 0 && bp_dy == 0) {
        return a.x < b.x;
    }

    if (ap_dy == 0) {
        return true;
    }

    if (bp_dy == 0) {
        return false;
    }

    double ap_inverse_slope = ap_dx / ap_dy;
    double bp_inverse_slope = bp_dx / bp_dy;

    return ap_inverse_slope > bp_inverse_slope;
}

// Normalize a double from the range [-1.0, 1.0] to range of 32-bit unsigned int.
//uint32_t norm_to_uint(double value) {
//    return static_cast<uint32_t>((value + 1) * 2147483647);
//}

// Normalize x and y to 32-bit unsigned integer range and compute Morton code.
//uint64_t morton_encode(const Point& p) {
//    return _pdep_u64(norm_to_uint(p.x), 0x5555555555555555)
//           | _pdep_u64(norm_to_uint(p.y), 0xaaaaaaaaaaaaaaaa);
//}

void argsort_x(const std::vector<Point>& points, std::vector<uint32_t>& index_array) {
    index_array.resize(points.size());
    std::iota(index_array.begin(), index_array.end(), 0);

    std::sort(
        index_array.begin(),
        index_array.end(),
        [&points](uint32_t i1, uint32_t i2) { return points[i1].x < points[i2].x; }
    );
}

void argsort_y(const std::vector<Point>& points, std::vector<uint32_t>& index_array) {
    index_array.resize(points.size());
    std::iota(index_array.begin(), index_array.end(), 0);

    std::sort(
        index_array.begin(),
        index_array.end(),
        [&points](uint32_t i1, uint32_t i2) { return points[i1].y < points[i2].y; }
    );
}

void compute_index_mask(const std::vector<uint32_t> indices, std::vector<bool>& mask, uint32_t n, bool init_mask = true) {
    if (init_mask) mask.resize(n);
    for (uint32_t i = 0; i < indices.size(); ++i) {
        mask[indices[i]] = true;
    }
}

void sort_points_and_indices(std::vector<Point>& points, std::vector<uint32_t>& point_indices,
                             std::vector<uint32_t>& permutation) {
    for (uint32_t i = 0; i < permutation.size(); ++i) {
        auto cur_index = i;
        while (i != permutation[cur_index]) {
            uint32_t next_index = permutation[cur_index];
            std::swap(points[cur_index], points[next_index]);
            std::swap(point_indices[cur_index], point_indices[next_index]);
            permutation[cur_index] = cur_index;
            cur_index = next_index;
        }
        permutation[cur_index] = cur_index;
    }
}

std::pair<uint32_t, uint32_t> build_2d_tree(const std::vector<Point>& points,
                                            std::vector<uint32_t>& point_indices,
                                            KDTree& tree, uint32_t depth) {
    std::vector<uint32_t> points_sort_idx;
    std::vector<Point> points_copy = points;

    if (depth % 2 == 0) {
        argsort_x(points_copy, points_sort_idx);
    } else {
        argsort_y(points_copy, points_sort_idx);
    }

    sort_points_and_indices(points_copy, point_indices, points_sort_idx);
    uint32_t mid = points_copy.size() / 2;
    uint32_t p = point_indices[mid];
    uint32_t max_depth = depth;

    if (mid > 0) {
        std::vector<uint32_t> point_indices_left;
        std::vector<Point> points_left;
        point_indices_left.reserve(mid);
        points_left.reserve(mid);

        for (uint32_t i = 0; i < mid; ++i) {
            point_indices_left.push_back(point_indices[i]);
        }

        for (uint32_t i = 0; i < mid; ++i) {
            points_left.push_back(points_copy[i]);
        }

        std::pair<uint32_t, uint32_t> left_result = build_2d_tree(points_left, point_indices_left, tree, depth + 1);
        max_depth = left_result.second;
        tree.nodes[p].set_left(left_result.first);
    } else {
        tree.nodes[p].set_left(p);
    }

    if (mid < points_copy.size() - 1) {
        std::vector<uint32_t> point_indices_right;
        std::vector<Point> points_right;
        uint32_t n_points_right = points_copy.size() - mid - 1;
        point_indices_right.reserve(n_points_right);
        points_right.reserve(n_points_right);

        for (uint32_t i = mid + 1; i < points_copy.size(); ++i) {
            point_indices_right.push_back(point_indices[i]);
        }

        for (uint32_t i = mid + 1; i < points_copy.size(); ++i) {
            points_right.push_back(points_copy[i]);
        }

        std::pair<uint32_t, uint32_t> right_result = build_2d_tree(points_right, point_indices_right, tree, depth + 1);
        if (right_result.second > max_depth) {
            max_depth = right_result.second;
        }
        tree.nodes[p].set_right(right_result.first);
    } else {
        tree.nodes[p].set_right(p);
    }

    return {p, max_depth};
}

KDTree::KDTree(const std::vector<Point>& points) {
    std::vector<uint32_t> point_indices(points.size());
    std::iota(point_indices.begin(), point_indices.end(), 0);
    nodes.resize(points.size());
    std::pair<uint32_t, uint32_t> result = build_2d_tree(points, point_indices, *this, 0);
    root = result.first;
    max_depth = result.second;
}

void compute_convex_hull_mask(const std::vector<Point>& points, std::vector<bool>& mask) {
    uint32_t lowest_p_index = 0;
    double lowest_y = points[0].y;
    for (uint32_t i = 1; i < points.size(); ++i) {
        if (points[i].y < lowest_y) {
            lowest_y = points[i].y;
            lowest_p_index = i;
        }
    }

    std::vector<uint32_t> index_array(points.size() - 1);
    uint32_t point_index = 0;
    while (point_index != lowest_p_index) {
        index_array[point_index] = point_index;
        ++point_index;
    }
    while (point_index != index_array.size()) {
        index_array[point_index] = point_index + 1;
        ++point_index;
    }

    std::sort(
        index_array.begin(),
        index_array.end(),
        [&points, lowest_p_index](uint32_t i1, uint32_t i2) {
            return compare_points_for_graham_scan(points[i1], points[i2], points[lowest_p_index]);
        }
    );

    std::vector<uint32_t> stack(points.size());
    stack[0] = lowest_p_index;
    uint32_t top = 0;

    for (uint32_t i : index_array) {
        while (top > 0 && points[i].orientation(points[stack[top - 1]], points[stack[top]]) < 0) {
            --top;
        }
        stack[++top] = i;
    }

    // truncate the stack to the actual size
    stack.resize(top + 1);

    compute_index_mask(stack, mask, points.size());
}

void gen_rand_points(std::vector<Point>& points, uint32_t n) {
    points.resize(n);
    for (uint32_t i = 0; i < n; ++i) {
        points[i] = {static_cast<double>(rand() * 2 - 1) / RAND_MAX,
                     static_cast<double>(rand() * 2 - 1) / RAND_MAX};
    }
}

void print_tree(const KDTree& tree) {
    for (uint32_t i = 0; i < tree.nodes.size(); ++i) {
        std::cout << i << ": " << tree.nodes[i].left_index << ", " << tree.nodes[i].right_index << std::endl;
    }
}

struct KDStackFrame {
    uint32_t point_idx;
    uint32_t depth;
    double diff;
    BBox bbox;

    KDStackFrame() : point_idx(0), depth(0), diff(0), bbox() {}

    void set(uint32_t point_index, uint32_t node_depth, double difference) {
        point_idx = point_index;
        depth = node_depth;
        diff = difference;
    }
};

struct KDStack {
    KDStack(uint32_t size) : stack(size), top(-1) {}

    void split_bbox_left(const BBox& bbox, double val, int axis, int frame_index) {
        stack[frame_index].bbox.set_left(bbox, val, axis);
    }

    void split_bbox_right(const BBox& bbox, double val, int axis, int frame_index) {
        stack[frame_index].bbox.set_right(bbox, val, axis);
    }

    void reset(uint32_t root) {
        top = 0;
        stack[0].point_idx = root;
        stack[0].depth = 0;
        stack[0].diff = HUGE_VAL;
        stack[0].bbox.reset();
    }

    const KDStackFrame& pop() {
        return stack[top--];
    }

    const bool empty() const {
        return top == -1;
    }

    void update(uint32_t idx_left, uint32_t idx_right,
                double diff, double val, int axis, bool use_diff=true) {
        double diff2 = use_diff ? diff : val;
        int prev_top = top + 1;
        uint32_t idx = stack[prev_top].point_idx;
        uint32_t next_depth = stack[prev_top].depth + 1;
        if (diff > 0) {
            if (idx_right != idx) {
                ++top;
                int frame_idx = (idx_left == idx) ? prev_top : prev_top + 1;
                split_bbox_right(stack[prev_top].bbox, val, axis, frame_idx);
                stack[frame_idx].set(idx_right, next_depth, HUGE_VAL);
            }
            if (idx_left != idx) {
                ++top;
                split_bbox_left(stack[prev_top].bbox, val, axis, prev_top);
                stack[prev_top].set(idx_left, next_depth, diff2);
            }
        } else {
            if (idx_left != idx) {
                ++top;
                int frame_idx = (idx_right == idx) ? prev_top : prev_top + 1;
                split_bbox_left(stack[prev_top].bbox, val, axis, frame_idx);
                stack[frame_idx].set(idx_left, next_depth, HUGE_VAL);
            }
            if (idx_right != idx) {
                ++top;
                split_bbox_right(stack[prev_top].bbox, val, axis, prev_top);
                stack[prev_top].set(idx_right, next_depth, diff2);
            }
        }
    }

    std::vector<KDStackFrame> stack;
    int top;
};

uint32_t find_vertex_approx(const KDTree& tree, const Point& query_point,
                            const std::vector<Point>& points, const Point& a, const Point& b,
                            const int orientation_sign, const Line& ab_line, KDStack& stack) {
    uint32_t best_point = 0;
    double best_dist = HUGE_VAL;
    stack.reset(tree.root);

    while (!stack.empty()) {
        const KDStackFrame& cur = stack.pop();

        if ((cur.diff != HUGE_VAL) && pow(cur.diff, 2) >= best_dist) {
            continue;
        }

        const double bbox_to_point_dist = cur.bbox.dist_to_point(query_point);
        if (bbox_to_point_dist >= best_dist) {
            continue;
        }

        const Point& p = points[cur.point_idx];

        if (p.orientation(a, b) * orientation_sign > 0) {
            double distance = p.distance_squared(query_point);
            if (distance < best_dist) {
                best_point = cur.point_idx;
                best_dist = distance;
            }
        } else if (bbox_to_point_dist != 0 && !cur.bbox.intersects(ab_line)) {
            continue;
        }

        const uint32_t axis = cur.depth % 2;
        const double diff = query_point[axis] - p[axis];

        stack.update(tree.nodes[cur.point_idx].left_index,
                     tree.nodes[cur.point_idx].right_index,
                     diff, p[axis], axis);
    }

    return best_point;
}

uint32_t find_vertex(const KDTree& tree, const std::vector<Point>& points,
                     const Point& a, const Point& b, const Point& c, KDStack& stack) {
    const int orientation_sign = c.orientation(a, b) > 0 ? -1 : 1;
    const PerpendicularUnitVector ab_perp(a, b, orientation_sign);
    const Line ab_line(a, b);
    const Point abc_center(b, c, ab_perp.line);
    const double abc_radius = a.distance(abc_center);
    const Point initial_query_point(abc_center.x + ab_perp.x * abc_radius,
                                    abc_center.y + ab_perp.y * abc_radius);
    uint32_t best_point = find_vertex_approx(tree, initial_query_point, points, a, b,
                                             orientation_sign, ab_line, stack);

    Point query_point(b, points[best_point], ab_perp.line);
    double best_dist = query_point.distance_squared(a);

    stack.reset(tree.root);

    while (!stack.empty()) {
        const KDStackFrame& cur = stack.pop();

        if (cur.diff != HUGE_VAL
            && pow(query_point[(cur.depth + 1) % 2] - cur.diff, 2) >= best_dist) {
            continue;
        }

        const double bbox_to_point_dist = cur.bbox.dist_to_point(query_point);
        if (bbox_to_point_dist >= best_dist) {
            continue;
        }

        const Point& p = points[cur.point_idx];

        if (p.orientation(a, b) * orientation_sign > 0) {
            double distance = p.distance_squared(query_point);
            if (distance < best_dist) {
                Point center(b, p, ab_perp.line);
                best_point = cur.point_idx;
                best_dist = center.distance_squared(p);
                query_point.x = center.x;
                query_point.y = center.y;
            }
        } else if (bbox_to_point_dist != 0 && !cur.bbox.intersects(ab_line)) {
            continue;
        }

        const uint32_t axis = cur.depth % 2;
        const double diff = query_point[axis] - p[axis];

        stack.update(tree.nodes[cur.point_idx].left_index,
                     tree.nodes[cur.point_idx].right_index,
                     diff, p[axis], axis, false);
    }

    return best_point;
}

struct Triangle {
    // vertex a (index into points vector)
    uint32_t a;

    // vertex b (index into points vector)
    uint32_t b;

    // vertex c (index into points vector)
    uint32_t c;

    // constructor
    Triangle(uint32_t v1, uint32_t v2, uint32_t v3) : a(v1), b(v2), c(v3) {}
    Triangle() : a(0), b(0), c(0) {}
};

Triangle find_initial_triangle(const std::vector<Point>& points, const std::vector<bool>& convex_hull_mask) {
    uint32_t v1 = 0;
    for (uint32_t i = 0; i < points.size(); ++i) {
        if (!convex_hull_mask[i]) {
            v1 = i;
            break;
        }
    }
    uint32_t v2 = 0;
    uint32_t v3 = 0;

    double min_distance = HUGE_VAL;

    for (uint32_t i = 0; i < points.size(); ++i) {
        if (i == v1) continue;
        double distance = points[v1].distance_squared(points[i]);
        if (distance < min_distance) {
            min_distance = distance;
            v2 = i;
        }
    }

    PerpendicularBisector v1_v2_perp(points[v1], points[v2]);
    double min_radius_squared = HUGE_VAL;

    for (uint32_t i = 0; i < points.size(); ++i) {
        if (i == v1 || i == v2) continue;

        Point center(points[v2], points[i], v1_v2_perp);
        double radius_squared = points[v1].distance_squared(center);
        if (radius_squared < min_radius_squared) {
            min_radius_squared = radius_squared;
            v3 = i;
        }
    }

    if (v1 < v2) {
        if (v1 < v3) {
            if (v2 < v3) {
                return Triangle(v1, v2, v3);
            }
            return Triangle(v1, v3, v2);
        }
        return Triangle(v3, v1, v2);
    }
    if (v1 < v3) {
        return Triangle(v2, v1, v3);
    }
    if (v2 < v3) {
        return Triangle(v2, v3, v1);
    }
    return Triangle(v3, v2, v1);
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (std::pair<T1, T2> const& pair) const {
        std::size_t h1 = std::hash<T1>()(pair.first);
        std::size_t h2 = std::hash<T2>()(pair.second);
        return h1 ^ h2;
    }
};

void triangulate(const std::vector<Point>& points, std::vector<Triangle>& triangles) {
    KDTree tree(points);
    KDStack stack(tree.max_depth + 1);
    std::vector<bool> convex_hull_mask;
    compute_convex_hull_mask(points, convex_hull_mask);

    Triangle initial_triangle = find_initial_triangle(points, convex_hull_mask);
    triangles.reserve(points.size() * 2 - 5);
    triangles.push_back(initial_triangle);

    std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t, pair_hash> edges_todo;

    // Insert initial edges with associated vertices
    edges_todo[{initial_triangle.a, initial_triangle.b}] = initial_triangle.c;
    if (!convex_hull_mask[initial_triangle.b] || !convex_hull_mask[initial_triangle.c]) {
        edges_todo[{initial_triangle.b, initial_triangle.c}] = initial_triangle.a;
    }
    edges_todo[{initial_triangle.a, initial_triangle.c}] = initial_triangle.b;

    //std::vector<double> times;
    //times.reserve(points.size() * 2 - 5);

    //auto t1 = std::chrono::high_resolution_clock::now();
    while (!edges_todo.empty()) {
        auto item = edges_todo.begin();
        std::pair<uint32_t, uint32_t> ab = item->first;
        uint32_t a, b;
        std::tie(a, b) = ab;
        uint32_t c = item->second;
        edges_todo.erase(item);

        //auto start = std::chrono::high_resolution_clock::now();
        uint32_t d = find_vertex(tree, points, points[a], points[b], points[c], stack);
        //auto end = std::chrono::high_resolution_clock::now();
        //times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        std::pair<uint32_t, uint32_t> ad = (a < d) ? std::make_pair(a, d) : std::make_pair(d, a);
        std::pair<uint32_t, uint32_t> bd = (b < d) ? std::make_pair(b, d) : std::make_pair(d, b);

        if ((edges_todo.erase(ad) == 0) && (!convex_hull_mask[a] || !convex_hull_mask[d])) {
            edges_todo[ad] = b;
        }

        if ((edges_todo.erase(bd) == 0) && (!convex_hull_mask[b] || !convex_hull_mask[d])) {
            edges_todo[bd] = a;
        }

        triangles.emplace_back(a, b, d);
    }
    //auto t2 = std::chrono::high_resolution_clock::now();
    //std::cout << "Triangulation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;
    //double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    //std::cout << "Average find_vertex time: " << avg_time << " ns" << std::endl;
}

// main function
int main(int argc, char** argv) {
    // parse command line arguments
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <number of points>" << std::endl;
        return EXIT_FAILURE;
    }

    uint32_t n_points = atoi(argv[1]);
    std::vector<Point> points;
    gen_rand_points(points, n_points);

    std::vector<Triangle> triangles;
    triangulate(points, triangles);

    return EXIT_SUCCESS;
}
