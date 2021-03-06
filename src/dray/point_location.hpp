// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef DRAY_POINT_LOCATION_HPP
#define DRAY_POINT_LOCATION_HPP

#include <dray/array.hpp>
#include <dray/linear_bvh_builder.hpp>
#include <dray/vec.hpp>

namespace dray
{

class PointLocator
{
  protected:
  BVH m_bvh;

  PointLocator ();

  public:
  PointLocator (BVH bvh);
  ~PointLocator ();

  struct Candidates
  {
    Array<int32> m_candidates;
    Array<int32> m_aabb_ids;
  };

  template <typename T>
  Candidates locate_candidates (const Array<Vec<T, 3>> points, int32 max_candidates);

  template <typename T>
  Candidates locate_candidates (const Array<Vec<T, 3>> points,
                                const Array<int32> active_idx,
                                int32 max_candidates);
};

} // namespace dray

#endif
