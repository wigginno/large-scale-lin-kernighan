// Values of sin(theta) for theta = pi/8, pi/4, 3pi/8, ..., 2pi
__constant float16 sin_theta = {
    0.38268343, 0.70710678, 0.92387953, 1.0, 0.92387953, 0.70710678, 0.38268343, 0.0,
    -0.38268343, -0.70710678, -0.92387953, -1.0, -0.92387953, -0.70710678, -0.38268343, 0.0
};

// Values of cos(theta) for theta = pi/8, pi/4, 3pi/8, ..., 2pi
__constant float16 cos_theta = {
    0.92387953, 0.70710678, 0.38268343, 0.0, -0.38268343, -0.70710678, -0.92387953, -1.0,
    -0.92387953, -0.70710678, -0.38268343, 0.0, 0.382683432, 0.70710678, 0.92387953, 1.0
};

// Interleave bits of x and y, so that x and y occupy even and odd bits, respectively
inline ulong16 interleave(ulong16 x, ulong16 y) {
    x = (x | (x << 16)) & 0x0000FFFF0000FFFFUL;
    x = (x | (x << 8))  & 0x00FF00FF00FF00FFUL;
    x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0FUL;
    x = (x | (x << 2))  & 0x3333333333333333UL;
    x = (x | (x << 1))  & 0x5555555555555555UL;

    y = (y | (y << 16)) & 0x0000FFFF0000FFFFUL;
    y = (y | (y << 8))  & 0x00FF00FF00FF00FFUL;
    y = (y | (y << 4))  & 0x0F0F0F0F0F0F0F0FUL;
    y = (y | (y << 2))  & 0x3333333333333333UL;
    y = (y | (y << 1))  & 0x5555555555555555UL;

    return x | (y << 1);
}

__kernel void compute_sample_points(__global float2* origins,
                                    __global uint* scaled_distances,
                                    __global ulong16* output) {
    int gid = get_global_id(0);

    // Scale origin to range of 32-bit integer
    // Convert origins to double2
    uint2 origin_scaled = convert_uint2_sat_rte((convert_double2(origins[gid]) + 1.0) * 2147483647);

    // Convert polar to cartesian coordinates and interleave x and y
    ulong16 x = clamp(convert_ulong16(round(cos_theta * scaled_distances[gid] + origin_scaled.x)), 0, 4294967295);
    ulong16 y = clamp(convert_ulong16(round(sin_theta * scaled_distances[gid] + origin_scaled.y)), 0, 4294967295);

    // Compute Morton codes
    ulong16 out = interleave(x, y);

    // Sort Morton codes to make things easier for CPU.
    ulong8 tmp0, tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;
    tmp0 = min(out.s13579BDF, out.s2468ACEF);
    tmp1 = max(out.s13579BDF, out.s2468ACEF);
    out.s13579BDF = tmp0;
    out.s2468ACEF = tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;
    tmp0 = min(out.s13579BDF, out.s2468ACEF);
    tmp1 = max(out.s13579BDF, out.s2468ACEF);
    out.s13579BDF = tmp0;
    out.s2468ACEF = tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;
    tmp0 = min(out.s13579BDF, out.s2468ACEF);
    tmp1 = max(out.s13579BDF, out.s2468ACEF);
    out.s13579BDF = tmp0;
    out.s2468ACEF = tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;
    tmp0 = min(out.s13579BDF, out.s2468ACEF);
    tmp1 = max(out.s13579BDF, out.s2468ACEF);
    out.s13579BDF = tmp0;
    out.s2468ACEF = tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;
    tmp0 = min(out.s13579BDF, out.s2468ACEF);
    tmp1 = max(out.s13579BDF, out.s2468ACEF);
    out.s13579BDF = tmp0;
    out.s2468ACEF = tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;
    tmp0 = min(out.s13579BDF, out.s2468ACEF);
    tmp1 = max(out.s13579BDF, out.s2468ACEF);
    out.s13579BDF = tmp0;
    out.s2468ACEF = tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;
    tmp0 = min(out.s13579BDF, out.s2468ACEF);
    tmp1 = max(out.s13579BDF, out.s2468ACEF);
    out.s13579BDF = tmp0;
    out.s2468ACEF = tmp1;
    tmp0 = min(out.even, out.odd);
    tmp1 = max(out.even, out.odd);
    out.even = tmp0;
    out.odd = tmp1;

    // Output sorted Morton codes
    output[gid] = out;
}
