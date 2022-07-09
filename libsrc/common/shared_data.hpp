#ifndef SHARED_DATA_HPP
#define SHARED_DATA_HPP

#include "config.hpp"
#include "mesh/mesh_block.hpp"
#include "utils/array_ref.hpp"

#include <vector>

namespace prl {

class shared_data_t;

////////////////////////////////////////////////////////////////////////////////
/// \brief Shared data class
///
/// You only need to add data here that you think other packages might want
/// to use.
////////////////////////////////////////////////////////////////////////////////

class shared_data_block_t {

  mesh_block_t & mesh_;

  int_t num_dims_ = 0;
  int_t num_cells_ = 0;
  int_t num_faces_ = 0;
  
public:

  //hydro
  std::vector<real_t> cell_density_;
  std::vector<real_t> cell_temperature_;
  std::vector<real_t> cell_cv_;

  shared_data_block_t(mesh_block_t & mesh) :
    mesh_(mesh)
  { reconfigure(); }

  void reconfigure();

  array_ref<real_t> resize_cell_density()
  {
    cell_density_.resize(num_cells_);
    return cell_density_;
  }
  
  array_ref<real_t> resize_cell_temperature()
  {
    cell_temperature_.resize(num_cells_);
    return cell_temperature_;
  }
  
  array_ref<real_t> resize_cv()
  {
    cell_cv_.resize(num_cells_);
    return cell_cv_;
  }
  
  friend class shared_data_t;
};

////////////////////////////////////////////////////////////////////////////////
/// \brief Shared data class
///
/// You only need to add data here that you think other packages might want
/// to use.
////////////////////////////////////////////////////////////////////////////////
class shared_data_t {
  std::vector<std::unique_ptr<shared_data_block_t>> blocks_;

public:

  void reserve(size_t n)
  { blocks_.reserve(n); }

  template<typename...Args>
  void emplace_back(Args && ... args)
  { 
    auto tmp = std::make_unique<shared_data_block_t>(std::forward<Args>(args)...);
    blocks_.emplace_back( std::move(tmp) );
  }

  size_t size() const
  { return blocks_.size(); }

  auto operator[](int_t i)
  { return blocks_[i].get(); }
  
  const auto operator[](int_t i) const
  { return blocks_[i].get(); }

  void reset(int_t i)
  { blocks_[i].reset(); }

  void post_amr();

};

} // namspace

#endif
