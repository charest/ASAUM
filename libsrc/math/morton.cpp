#include "morton.hpp"

namespace prl {

std::size_t morton_number(int depth, unsigned x, unsigned y)
{
  std::size_t mask = 1 << (depth - 1);
  std::size_t result = 0;
  for (int b = depth; b--; ) {
    result |= y & mask; result <<= 1;
    result |= x & mask; mask >>= 1;
  }
  return result;
}

std::size_t morton_number(int depth, unsigned x, unsigned y, unsigned z)
{
  std::size_t mask = 1 << (depth - 1);
  std::size_t result = 0;
  for (int b = depth; b--; ) {
    result |= x & mask; result <<= 1;
    result |= y & mask; result <<= 1;
    result |= z & mask; mask >>= 1;
  }
  return result;
}

} // namespace
