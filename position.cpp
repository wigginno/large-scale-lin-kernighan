#include <immintrin.h>

// Normalize a 32-bit unsigned integer to a floating point number in the range [-1.0, 1.0).
inline float norm(uint32_t value) {
    return static_cast<double>(value) / 2147483647 - 1.0;
}

// Denormalize a floating point number from the range [-1.0, 1.0) back to a 32-bit unsigned integer.
inline uint32_t denorm(float value) {
    return (static_cast<double>(value) + 1.0) * 2147483647;
}

// Interleave bits of x and y, so that x and y occupy even and odd bits, respectively.
inline uint64_t interleave(uint32_t x, uint32_t y) {
    return _pdep_u64(x, 0x5555555555555555) | _pdep_u64(y, 0xaaaaaaaaaaaaaaaa);
}

// Place even bits of z in x and odd bits in y
inline void deinterleave(uint64_t z, uint32_t &x, uint32_t &y) {
    x = _pext_u64(z, 0x5555555555555555);
    y = _pext_u64(z, 0xaaaaaaaaaaaaaaaa);
}

inline void morton_to_position(uint64_t z, float &x, float &y) {
    uint32_t x_scaled, y_scaled;
    deinterleave(z, x_scaled, y_scaled);
    x = norm(x_scaled);
    y = norm(y_scaled);
}

inline uint64_t position_to_morton(float x, float y) {
    uint32_t x_scaled = denorm(x);
    uint32_t y_scaled = denorm(y);
    return interleave(x_scaled, y_scaled);
}
