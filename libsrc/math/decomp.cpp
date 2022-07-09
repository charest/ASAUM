#include "decomp.hpp"

using std::size_t;

namespace prl {

std::vector<size_t> sieve(size_t n)
{
  std::vector<bool> prime(n+1, true);
  std::vector<size_t> ret;

  for (size_t p = 2; p*p <= n; p++) {
    if (prime[p]) {
      for (size_t i = 2*p; i <= n; i += p) {
        prime[i] = false;
      }
    }
  }

  for (size_t p = 2; p <= n; p++) {
    if (prime[p])
      ret.push_back(p);
  }

  return ret;
}

} // namespace
