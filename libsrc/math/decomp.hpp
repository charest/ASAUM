#ifndef DECOMP_HPP
#define DECOMP_HPP

#include <array>
#include <algorithm>
#include <vector>

namespace prl {

std::vector<std::size_t> sieve(std::size_t n);

template<typename T>
std::vector<T> factor(size_t np)
{
  std::vector<T> facs;
  auto primes = sieve(np);

  size_t p = 0;
  for (size_t i = 0; i < primes.size(); i++) {
    if (np % primes[p] == 0) {
      facs.push_back(primes[p]);
      np = np / primes[p];
    } else {
      p++;
    }
  }

  std::sort(facs.begin(), facs.end());
  std::reverse(facs.begin(), facs.end());

  return facs;
}


template<typename T>
std::vector<T> decomp(std::vector<T> n, std::size_t np)
{
  auto d = n.size();
  if (d==0) return {};

  std::vector<T> parts(d, 1);

  auto facs = factor<T>(np);

  // greedy decomposition
  for (auto fac : facs) {
    auto maxind = std::distance(n.begin(), std::max_element(n.begin(), n.end()));
    parts[maxind] *= fac;
    n[maxind] /= fac;
  }

  return parts;
}

} // namespace

#endif
