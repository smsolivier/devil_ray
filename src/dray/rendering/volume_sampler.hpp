// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef DRAY_VOLUME_SAMPLER_HPP
#define DRAY_VOLUME_SAMPLER_HPP

#include <dray/bounds_map.hpp>
#include <dray/collection.hpp>
#include <dray/color_map.hpp>
#include <dray/ray.hpp>
#include <dray/ray_hit.hpp>
#include <dray/rendering/domain_sampler.hpp>
#include <dray/rendering/point_light.hpp>

namespace dray
{

class VolumeSampler
{
protected:
  ColorMap m_color_map;
  Collection m_collection;
  std::string m_field;
  Range m_field_range;
  std::map<int32, DomainSampler> m_domains;
  BoundsMap m_bounds_map;

public:
  VolumeSampler() = delete;
  VolumeSampler(Collection &collection);
  ~VolumeSampler();

  void sample(Array<Ray> &rays);

  /// set the input data collection
  void input(Collection &collection);

  void field(const std::string field);
  std::string field() const;

protected:
  void init();
};


} // namespace dray
#endif
