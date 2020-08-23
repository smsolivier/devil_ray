#ifndef DRAY_BOUNDS_MAP_HPP
#define DRAY_BOUNDS_MAP_HPP

#include <dray/aabb.hpp>
#include <dray/bvh.hpp>
#include <map>

namespace dray
{

class BoundsMap
{
public:
  void clear();

  void add_block(int32 id, const AABB<3> &bounds);

  // returns -1 if not found
  int get_rank(const int32 &block_id);

  void build();

protected:
  std::map<int32, AABB<3>> m_bounds; // map<domain_id, bounds>
  std::map<int32, int32> m_rank_map;   // map<domain_id, rank>
  BVH m_bvh;
};

} // namespace dray
#endif
