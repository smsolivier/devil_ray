// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include <dray/rendering/volume_sampler.hpp>
#include <dray/rendering/device_framebuffer.hpp>
#include <dray/rendering/colors.hpp>
#include <dray/rendering/surface.hpp>
#include <dray/rendering/volume_shader.hpp>

#include <dray/filters/mesh_boundary.hpp>

#include <dray/dispatcher.hpp>
#include <dray/array_utils.hpp>
#include <dray/error_check.hpp>
#include <dray/device_color_map.hpp>

#include <dray/utils/data_logger.hpp>
#include <dray/utils/timer.hpp>

#include <dray/GridFunction/device_mesh.hpp>
#include <dray/GridFunction/device_field.hpp>

namespace dray
{

namespace detail
{

struct RaySegmentFunctor
{
  Array<Ray> *m_rays;
  Array<Vec<Float,2>> m_start_distance;
  RaySegmentFunctor(Array<Ray> *rays)
    : m_rays(rays)
  {
  }

  template<typename TopologyType>
  void operator()(TopologyType &topo)
  {

  }

};

void write_segments(Array<Vec<Float,2>> &segments, Array<Ray> &rays)
{
  Framebuffer fb(1000,1000);
  fb.clear();
  const int size = rays.size();

  Vec<Float,2> *seg_ptr = segments.get_host_ptr();
  Ray *ray_ptr = rays.get_host_ptr();
  float32 * depths_ptr = fb.depths().get_host_ptr();

  for(int i = 0 ;  i < size; ++i)
  {
    float32 depth = 0.f;
    Vec<Float,2> seg = seg_ptr[i];
    if(seg[0] >= 0.f)
    {
      depth = seg[1] - seg[0];
      //std::cout<<"Depths = "<<depth<<" "<<seg<<"\n";
    }
    depths_ptr[ray_ptr[i].m_pixel_id] = depth;
  }

  fb.save_depth("segment");

}

Array<Vec<Float,3>> origin_to_points(Array<Ray> &rays)
{
  const int32 size = rays.size();
  const Ray *ray_ptr = rays.get_device_ptr_const();

  Array<Vec<Float,3>> points;
  points.resize(size);

  Vec<Float,3> *points_ptr = points.get_device_ptr();

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, size), [=] DRAY_LAMBDA (int32 i)
  {
    points_ptr[i] = ray_ptr[i].m_orig;
  });

  DRAY_ERROR_CHECK();
  return points;
}

void bump_min_dist(Array<Ray> &rays, Array<RayHit> &hits, const float32 eps)
{
  const int32 size = rays.size();
  Ray *ray_ptr = rays.get_device_ptr();
  const RayHit *hit_ptr = hits.get_device_ptr_const();

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, size), [=] DRAY_LAMBDA (int32 i)
  {
    Ray ray = ray_ptr[i];
    RayHit hit = hit_ptr[i];
    if(hit.m_hit_idx != -1)
    {
      ray.m_near = hit.m_dist + eps;
    }
    ray_ptr[i] = ray;
  });

  DRAY_ERROR_CHECK();
}

Array<Vec<Float,2>> calc_segment(Array<RayHit> &first_hits,
                                 Array<RayHit> &second_hits,
                                 Array<Location> &locs)
{
  const int32 size = first_hits.size();
  const RayHit *fhit_ptr = first_hits.get_device_ptr_const();
  const RayHit *shit_ptr = second_hits.get_device_ptr_const();
  const Location *loc_ptr = locs.get_device_ptr_const();

  Array<Vec<Float,2>> segments;
  segments.resize(size);
  Vec<Float,2> *seg_ptr = segments.get_device_ptr();


  RAJA::forall<for_policy>(RAJA::RangeSegment(0, size), [=] DRAY_LAMBDA (int32 i)
  {
    Vec<Float,2> res;
    res[0] = -1.f;
    res[1] = -1.f;

    Location loc = loc_ptr[i];
    RayHit first = fhit_ptr[i];
    RayHit second = shit_ptr[i];

    bool valid = true;
    if(loc.m_cell_id != -1)
    {
      // The ray origin is inside the mesh
      // so the segment begins at the origin of
      // the ray
      //std::cout<<"Inside "<<"\n";
      res[0] = 0.f;
      if(first.m_hit_idx != -1)
      {
        res[1] = first.m_dist;
      }
      else
      {
        valid = false;
      }
    }
    else if(first.m_hit_idx != -1 && second.m_hit_idx != -1)
    {
      // so we were not inside the mesh, so the next two
      // intersections should be what we are looking for
      res[0] = first.m_dist;
      res[1] = second.m_dist;
      //std::cout<<"Boom "<<first<<" "<<second<<"\n";
    }
    else
    {
      valid = false;
    }

    if(!valid)
    {
      res[0] = -1;
      //std::cout<<"V";
    }

    seg_ptr[i] = res;

  });

  DRAY_ERROR_CHECK();
  return segments;
}

} // namespace detail

// ------------------------------------------------------------------------
VolumeSampler::VolumeSampler(DataSet &data_set)
  : m_data_set(data_set)
{
  init();
}

// ------------------------------------------------------------------------
VolumeSampler::~VolumeSampler()
{
}

// ------------------------------------------------------------------------
void
VolumeSampler::input(DataSet &data_set)
{
  m_data_set = data_set;
  init();
}

// ------------------------------------------------------------------------

void
VolumeSampler::field(const std::string field)
{
  m_field = field;
  std::vector<Range> ranges = m_data_set.field(m_field)->range();
  if(ranges.size() != 0)
  {
    DRAY_ERROR("Volume Sampler field must be a scalar");
  }
  m_field_range = ranges[0];
}

std::string
VolumeSampler::field() const
{
  return m_field;
}

void VolumeSampler::init()
{
  MeshBoundary boundary;
  m_boundary = boundary.execute(m_data_set);
}

void VolumeSampler::sample(Array<Ray> &rays)
{

  Collection col;
  col.add_domain(m_data_set);

  Collection boundary_col;
  boundary_col.add_domain(m_boundary);

  Surface surface(boundary_col);
  surface.active_domain(0);
  Array<RayHit> first_hits = surface.nearest_hit(rays);


  TopologyBase *topo = m_data_set.topology();
  Array<Vec<Float,3>> origins = detail::origin_to_points(rays);
  Array<Location> locs = topo->locate(origins);

  AABB<3> bounds = topo->bounds();
  std::cout<<" BOUnds "<<bounds<<"\n";
  std::cout<<" mag "<<(bounds.max() - bounds.min()).magnitude()<<"\n";
  const float32 eps = (bounds.max() - bounds.min()).magnitude()*1e-4;
  std::cout<<"Eps "<<eps<<"\n";

  Array<Ray> ray_copy;
  array_copy(ray_copy, rays);

  detail::bump_min_dist(ray_copy, first_hits,  eps);
  Array<RayHit> second_hits = surface.nearest_hit(ray_copy);


  Array<Vec<Float,2>> segments;
  segments = detail::calc_segment(first_hits, second_hits, locs);

  detail::write_segments(segments, rays);
  //topo->locate(
  //TopologyBase *topo = m_boundary.topology();
  //detail::RaySegmentFunctor func(&rays);
  //dispatch_2d(topo, func);
}

// ------------------------------------------------------------------------
#if 0
Array<VolumePartial>
Volume::integrate(Array<Ray> &rays, Array<PointLight> &lights)
{

  DataSet data_set = m_collection.domain(m_active_domain);
  if(m_field == "")
  {
    DRAY_ERROR("Field never set");
  }
  if(!m_color_map.range_set())
  {
    m_color_map.scalar_range(m_field_range);
  }

  TopologyBase *topo = data_set.topology();
  FieldBase *field = data_set.field(m_field);

  detail::IntegratePartialsFunctor func(&rays,
                                        lights,
                                        m_color_map,
                                        m_samples,
                                        m_bounds,
                                        m_use_lighting);
  dispatch_3d(topo, field, func);
  return func.m_partials;
}
// ------------------------------------------------------------------------

void Volume::samples(int32 num_samples)
{
  m_samples = num_samples;
}

// ------------------------------------------------------------------------

void Volume::use_lighting(bool do_it)
{
  m_use_lighting = do_it;
}

#endif
// ------------------------------------------------------------------------
} // namespace dray
