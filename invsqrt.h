#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define USE_OWN_INTRINSICS 0

#if USE_OWN_INTRINSICS
__forceinline int __clz (uint32_t a)
{
    int r;
    __asm__ ("bsrl %1,%0\n\t" : "=r"(r): "r"(a));
    return 31 - r;
}

uint32_t __umulhi (uint32_t a, uint32_t b)
{
    uint32_t r;
    __asm__ ("movl %1,%%eax\n\tmull %2\n\tmovl %%edx,%0\n\t"
             : "=r"(r) : "r"(a), "r"(b) : "eax", "edx");
    return r;
}
#else // USE_OWN_INTRINSICS
int __clz (uint32_t a)
{
    uint32_t r = 32;
    if (a >= 0x00010000) { a >>= 16; r -= 16; }
    if (a >= 0x00000100) { a >>=  8; r -=  8; }
    if (a >= 0x00000010) { a >>=  4; r -=  4; }
    if (a >= 0x00000004) { a >>=  2; r -=  2; }
    r -= a - (a & (a >> 1));
    return r;
}

uint32_t __umulhi (uint32_t a, uint32_t b)
{
    return (uint32_t)(((uint64_t)a * b) >> 32);
}
#endif // USE_OWN_INTRINSICS

/*
 * For each sub-interval in [1, 4), use an 8-bit approximation r to reciprocal
 * square root. To speed up subsequent Newton-Raphson iterations, each entry in
 * the table combines two pieces of information: The least-significant 10 bits
 * store 3*r, the most-significant 22 bits store r**3, rounded from 24 down to
 * 22 bits such that accuracy is optimized.
 */
uint32_t rsqrt_tab [96] = 
{
    0xfa0bdefa, 0xee6af6ee, 0xe5effae5, 0xdaf27ad9,
    0xd2eff6d0, 0xc890aec4, 0xc10366bb, 0xb9a71ab2,
    0xb4da2eac, 0xadce7ea3, 0xa6f2b29a, 0xa279a694,
    0x9beb568b, 0x97a5c685, 0x9163067c, 0x8d4fd276,
    0x89501e70, 0x8563da6a, 0x818ac664, 0x7dc4fe5e,
    0x7a122258, 0x7671be52, 0x72e44a4c, 0x6f68fa46,
    0x6db22a43, 0x6a52623d, 0x67041a37, 0x65639634,
    0x622ffe2e, 0x609cba2b, 0x5d837e25, 0x5bfcfe22,
    0x58fd461c, 0x57838619, 0x560e1216, 0x53300a10,
    0x51c72e0d, 0x50621a0a, 0x4da48204, 0x4c4c2e01,
    0x4af789fe, 0x49a689fb, 0x485a11f8, 0x4710f9f5,
    0x45cc2df2, 0x448b4def, 0x421505e9, 0x40df5de6,
    0x3fadc5e3, 0x3e7fe1e0, 0x3d55c9dd, 0x3d55d9dd,
    0x3c2f41da, 0x39edd9d4, 0x39edc1d4, 0x38d281d1,
    0x37bae1ce, 0x36a6c1cb, 0x3595d5c8, 0x3488f1c5,
    0x3488fdc5, 0x337fbdc2, 0x3279ddbf, 0x317749bc,
    0x307831b9, 0x307879b9, 0x2f7d01b6, 0x2e84ddb3,
    0x2d9005b0, 0x2d9015b0, 0x2c9ec1ad, 0x2bb0a1aa,
    0x2bb0f5aa, 0x2ac615a7, 0x29ded1a4, 0x29dec9a4,
    0x28fabda1, 0x2819e99e, 0x2819ed9e, 0x273c3d9b,
    0x273c359b, 0x2661dd98, 0x258ad195, 0x258af195,
    0x24b71192, 0x24b6b192, 0x23e6058f, 0x2318118c,
    0x2318718c, 0x224da189, 0x224dd989, 0x21860d86,
    0x21862586, 0x20c19183, 0x20c1b183, 0x20001580
};

/* This function computes the reciprocal square root of its 16.16 fixed-point 
 * argument. After normalization of the argument if uses the most significant
 * bits of the argument for a table lookup to obtain an initial approximation 
 * accurate to 8 bits. This is followed by two Newton-Raphson iterations with
 * quadratic convergence. Finally, the result is denormalized and some simple
 * rounding is applied to maximize accuracy.
 *
 * To speed up the first NR iteration, for the initial 8-bit approximation r0
 * the lookup table supplies 3*r0 along with r0**3. A first iteration computes
 * a refined estimate r1 = 1.5 * r0 - x * r0**3. The second iteration computes
 * the final result as r2 = 0.5 * r1 * (3 - r1 * (r1 * x)).
 *
 * The accuracy for all arguments in [0x00000001, 0xffffffff] is as follows: 
 * 639 results are too small by one ulp, 1457 results are too big by one ulp.
 * A total of 2096 results deviate from the correctly rounded result.
 */
uint32_t fxrsqrt (uint32_t a)
{
    uint32_t s, r, t, scal;

    /* handle special case of zero input */
    if (a == 0) return ~a;
    /* normalize argument */
    scal = __clz (a) & 0xfffffffe;
    a = a << scal;
    /* initial approximation */
    t = rsqrt_tab [(a >> 25) - 32];
    /* first NR iteration */
    r = (t << 22) - __umulhi (t, a);
    /* second NR iteration */
    s = __umulhi (r, a);
    s = 0x30000000 - __umulhi (r, s);
    r = __umulhi (r, s);
    /* denormalize and round result */
    r = ((r >> (18 - (scal >> 1))) + 1) >> 1;
    return r;
}

/* reference implementation, 16.16 reciprocal square root of non-zero argment */
uint32_t ref_fxrsqrt (uint32_t a)
{
    double arg = a / 65536.0;
    double rsq = sqrt (1.0 / arg);
    uint32_t r = (uint32_t)(rsq * 65536.0 + 0.5);
    return r;
}
