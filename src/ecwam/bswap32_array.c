#include <stdint.h>
#include <stddef.h>

/* Byte-reverse every 4-byte element in buf[0..n-1] in-place.
 * __builtin_bswap32 maps to a single BSWAP instruction; the loop
 * auto-vectorises to VPSHUFB on AVX2 targets and is parallelised
 * with OpenMP across available threads. */
void bswap32_array(float *buf, int n)
{
    uint32_t *p = (uint32_t *)buf;
    size_t    i, len = (size_t)n;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < len; i++)
        p[i] = __builtin_bswap32(p[i]);
}
