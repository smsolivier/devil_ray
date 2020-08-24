// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef DRAY_DOMAIN_SAMPLER_HPP
#define DRAY_DOMAIN_SAMPLER_HPP

#include <dray/collection.hpp>
#include <dray/color_map.hpp>
#include <dray/ray.hpp>
#include <dray/ray_hit.hpp>
#include <dray/rendering/volume_partial.hpp>
#include <dray/rendering/point_light.hpp>

namespace dray
{

class DomainSampler
{
protected:
  ColorMap m_color_map;
  DataSet m_data_set;
  DataSet m_boundary;
  std::string m_field;
  Range m_field_range;

public:
  DomainSampler();
  DomainSampler(const DataSet &data_set);
  ~DomainSampler();

  //Array<VolumePartial> integrate(Array<Ray> &rays, Array<PointLight> &lights);
  void sample(Array<Ray> &rays);

  /// set the input data set
  void input(DataSet &data_set);

  void field(const std::string field);
  void range(const Range &range);
  std::string field() const;

protected:
  void init();
};


} // namespace dray
#endif
