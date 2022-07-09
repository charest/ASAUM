#ifndef MORTON_HPP
#define MORTON_HPP

#include <cstddef>

namespace prl {

std::size_t morton_number(int depth, unsigned x, unsigned y);
std::size_t morton_number(int depth, unsigned x, unsigned y, unsigned z);

} // namespace

#endif
