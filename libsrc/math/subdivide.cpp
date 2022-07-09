#include <stddef.h>
#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// is an id owned by this rank
////////////////////////////////////////////////////////////////////////////////
bool owned_by( const std::vector<size_t> & dist, size_t id, size_t owner)
{
  if (id < dist[owner]) return false;
  else if (id < dist[owner+1]) return true;
  else return false;
}

} // namespace
