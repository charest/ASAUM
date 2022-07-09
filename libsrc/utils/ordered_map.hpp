#ifndef ORDERED_MAP_HPP
#define ORDERED_MAP_HPP

namespace prl {

template<typename Key, typename T>
struct ordered_map {
  
  using map_type = std::map<Key, T>;
  
  map_type map;
  std::vector<typename map_type::iterator> order;

  auto insert(const std::pair<Key,T> & val)
  {
    auto res = map.insert(val);
    if (res.second) order.emplace_back(res.first);
    return res;
  }

  template<typename...Args>
  auto emplace(Args &&... args)
  {
    auto res = map.emplace(std::forward<Args>(args)...);
    if (res.second) order.emplace_back(res.first);
    return res;
  }

  T & operator[](const Key & key) {
    auto it = map.find(key);
    if (it == map.end()) {
      T newval;
      auto res = insert({key, newval});
      return res.first.second;
    }
    else {
      return it->second;
    }
  }

  bool empty() const { return map.empty(); }

  size_t size() const { return map.size(); }

  void clear() {
    map.clear();
    order.clear();
  }
};

} // namespace

#endif
