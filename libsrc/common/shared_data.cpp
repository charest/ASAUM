#include "shared_data.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Reconfigure mesh
////////////////////////////////////////////////////////////////////////////////
void shared_data_block_t::reconfigure()
{
  num_dims_ = mesh_.num_dims();
  num_cells_ = mesh_.num_all_cells();
  num_faces_ = mesh_.num_all_faces();
}

////////////////////////////////////////////////////////////////////////////////
/// Post amr
////////////////////////////////////////////////////////////////////////////////
void shared_data_t::post_amr()
{
  auto last = std::remove_if(
      blocks_.begin(),
      blocks_.end(),
      [](const auto & a) { return !a; }
  );
  blocks_.erase(last, blocks_.end());

  std::sort(
      blocks_.begin(),
      blocks_.end(),
      [](const auto & a, const auto & b)
      { return a->mesh_.id() < b->mesh_.id(); }
  );

  for (auto & b : blocks_)
    b->reconfigure();
}

} // namespace
