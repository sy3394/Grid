/*
 * Diagnostic for crc32 symbol shadowing in the Aurora oneAPI/SYCL link.
 *
 * The standard CRC-32 of the ASCII string "123456789" is 0xcbf43926.
 * If this program prints anything else, the crc32 symbol in the link is
 * not zlib's, and every SciDAC checksum formed through it is wrong.
 *
 * Build with the same environment/flags/link line as the tracker binary
 * (after sourcing sourceme.sh), e.g.:
 *
 *   icpx $CXXFLAGS crc32-diag.cc $LDFLAGS -lz -o crc32-diag
 *
 * then:
 *
 *   ./crc32-diag
 *   ldd crc32-diag | grep -i libz
 *   for lib in $(ldd crc32-diag | awk '/=>/ {print $3}'); do
 *     nm -D --defined-only "$lib" 2>/dev/null | grep -qw crc32 && echo "$lib exports crc32";
 *   done
 */
#include <zlib.h>
#include <cstdio>

int main(void)
{
  const unsigned char check[] = "123456789";
  unsigned long c = crc32(0L, check, 9);
  printf("crc32(\"123456789\") = %08lx (expect cbf43926)\n", c);
  if (c == 0xcbf43926UL) {
    printf("OK: crc32 resolves to a standard zlib implementation\n");
    return 0;
  }
  printf("FAIL: crc32 symbol is shadowed by a non-zlib implementation\n");
  return 1;
}
