// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <dray/rendering/scalar_renderer.hpp>

#include <dray/dray.hpp>
#include <dray/error.hpp>
#include <dray/error_check.hpp>
#include <dray/rendering/volume.hpp>
#include <dray/utils/data_logger.hpp>
#include <dray/policies.hpp>

#include <memory>
#include <vector>

#include <apcomp/scalar_compositor.hpp>

#ifdef DRAY_MPI_ENABLED
#include <mpi.h>
#endif

namespace dray
{

namespace
{

void calc_offsets(Collection &collection, std::vector<int32> &offsets)
{
  const int size = collection.size();
  offsets.resize(size);
  int total_cells = 0;
  for(int i = 0; i < size; ++i)
  {
    offsets[i] = total_cells;

    total_cells += collection.domain(i).topology()->cells();
  }

#ifdef DRAY_MPI_ENABLED
  MPI_Comm comm = MPI_Comm_f2c(dray::mpi_comm());
  int *global_counts = new int[dray::mpi_size()];
  MPI_Allgather(&total_cells, 1, MPI_INT,global_counts, 1, MPI_INT, comm);
  int rank = dray::mpi_rank();
  int global_offset = 0;

  if(global_counts[dray::mpi_size() - 1] < 0)
  {
    DRAY_ERROR("Zone count overflow: blame matt");
  }

  for(int i = 0; i < rank; ++i)
  {
    global_offset += global_counts[i];
  }

  delete[] global_counts;

  for(int i = 0; i < size; ++i)
  {
    offsets[i] += global_offset;
  }
#endif
}

ScalarBuffer convert(apcomp::ScalarImage &image, std::vector<std::string> &names)
{
  const int num_fields = names.size();

  const int dx  = image.m_bounds.m_max_x - image.m_bounds.m_min_x + 1;
  const int dy  = image.m_bounds.m_max_y - image.m_bounds.m_min_y + 1;
  const int size = dx * dy;

  ScalarBuffer result(dx, dy, nan<float32>());

  std::vector<float*> buffers;
  for(int i = 0; i < num_fields; ++i)
  {
    result.add_field(names[i]);
    float* buffer = result.m_scalars[names[i]].get_host_ptr();
    buffers.push_back(buffer);
  }

  int32 *zone_id_ptr = result.m_zone_ids.get_host_ptr();

  // this is fields + zone_id for each pixel
  const unsigned char *loads = &image.m_payloads[0];
  const size_t payload_size = image.m_payload_bytes;

#ifdef DRAY_OPENMP_ENABLED
    #pragma omp parallel for
#endif
  for(size_t x = 0; x < size; ++x)
  {
    for(int i = 0; i < num_fields; ++i)
    {
      const size_t offset = x * payload_size + i * sizeof(float);
      memcpy(&buffers[i][x], loads + offset, sizeof(float));
    }

    const size_t zone_offset = x * payload_size + num_fields * sizeof(float);
    memcpy(zone_id_ptr+x, loads + zone_offset, sizeof(int32));
  }

  result.m_depths.resize(size);
  float* dbuffer = result.m_depths.get_host_ptr();
  memcpy(dbuffer, &image.m_depths[0], sizeof(float) * size);

  return result;
}

apcomp::ScalarImage *
convert(ScalarBuffer &result, const std::vector<std::string> &names)
{
  const int num_fields = names.size();
  // payload is fields + zone id
  const int payload_size = num_fields * sizeof(float) + sizeof(int32);

  apcomp::Bounds bounds;
  bounds.m_min_x = 1;
  bounds.m_min_y = 1;
  bounds.m_max_x = result.m_width;
  bounds.m_max_y = result.m_height;
  const size_t size = result.m_width * result.m_height;

  apcomp::ScalarImage *image = new apcomp::ScalarImage(bounds, payload_size);
  unsigned char *loads = &image->m_payloads[0];

  const float* dbuffer = result.m_depths.get_host_ptr();
  memcpy(&image->m_depths[0], dbuffer, sizeof(float) * size);

  // copy scalars into payload
  std::vector<float*> buffers;
  for(auto field : names)
  {
    float* buffer = result.m_scalars[field].get_host_ptr();
    buffers.push_back(buffer);
  }

  const int32 *zone_id_ptr = result.m_zone_ids.get_host_ptr_const();

#ifdef DRAY_OPENMP_ENABLED
    #pragma omp parallel for
#endif
  for(size_t x = 0; x < size; ++x)
  {
    for(int i = 0; i < num_fields; ++i)
    {
      const size_t offset = x * payload_size + i * sizeof(float);
      memcpy(loads + offset, &buffers[i][x], sizeof(float));
    }
    const size_t zone_id_offset = x * payload_size + num_fields * sizeof(float);
    memcpy(loads + zone_id_offset, zone_id_ptr+x, sizeof(int32));
  }
  return image;
}

} // namespace


ScalarRenderer::ScalarRenderer()
  : m_traceable(nullptr)
{
}

ScalarRenderer::ScalarRenderer(std::shared_ptr<Traceable> traceable)
{
  set(traceable);
}

void ScalarRenderer::set(std::shared_ptr<Traceable> traceable)
{
  m_traceable = traceable;
  calc_offsets(m_traceable->collection(), m_offsets);
}

void ScalarRenderer::field_names(const std::vector<std::string> &field_names)
{
  m_field_names = field_names;
}

ScalarBuffer
ScalarRenderer::render(Camera &camera)
{
  if(m_traceable == nullptr)
  {
    DRAY_ERROR("ScalarRenderer: traceable never set");
  }

  Array<Ray> rays;
  camera.create_rays (rays);

  ScalarBuffer scalar_buffer(camera.get_width(),
                             camera.get_height(),
                             nan<float32>());

  const int32 buffer_size = camera.get_width() * camera.get_height();
  const int32 field_size = m_field_names.size();

  scalar_buffer.m_depths.resize(buffer_size);

  const int domains = m_traceable->num_domains();

  for(int d = 0; d < domains; ++d)
  {
    DRAY_INFO("Tracing scalar domain "<<d);
    m_traceable->active_domain(d);
    const int domain_offset = m_offsets[d];

    Array<RayHit> hits = m_traceable->nearest_hit(rays);
    float32 *depth_ptr = scalar_buffer.m_depths.get_device_ptr();
    int32 *zone_id_ptr = scalar_buffer.m_zone_ids.get_device_ptr();

    const Ray *ray_ptr = rays.get_device_ptr_const ();
    const RayHit *hit_ptr = hits.get_device_ptr_const ();

    RAJA::forall<for_policy> (RAJA::RangeSegment (0, hits.size ()), [=] DRAY_LAMBDA (int32 ii)
    {
      const RayHit &hit = hit_ptr[ii];
      const Ray &ray = ray_ptr[ii];
      if (hit.m_hit_idx > -1)
      {
        const int32 pid = ray.m_pixel_id;
        depth_ptr[pid] = hit.m_dist;
        zone_id_ptr[pid] = hit.m_hit_idx + domain_offset;
      }
    });

    for(int32 i = 0; i < field_size; ++i)
    {
      std::string field = m_field_names[i];
      m_traceable->field(field);
      Array<Fragment> fragments = m_traceable->fragments(hits);
      if(!scalar_buffer.has_field(field))
      {
        scalar_buffer.add_field(field);
      }
      Array<float32> buffer = scalar_buffer.m_scalars[field];;
      float32 *buffer_ptr = buffer.get_device_ptr();

      const Fragment *frag_ptr = fragments.get_device_ptr_const ();
      RAJA::forall<for_policy> (RAJA::RangeSegment (0, hits.size ()), [=] DRAY_LAMBDA (int32 ii)
      {
        const RayHit &hit = hit_ptr[ii];
        const Fragment &frag = frag_ptr[ii];
        const Ray &ray = ray_ptr[ii];

        if (hit.m_hit_idx > -1)
        {
          const int32 pid = ray.m_pixel_id;
          buffer_ptr[pid] = frag.m_scalar;
        }

      });

    }

    ray_max(rays, hits);
  }


#ifdef DRAY_MPI_ENABLED
  apcomp::PayloadCompositor compositor;
  apcomp::ScalarImage *pimage = convert(scalar_buffer, m_field_names);
  compositor.AddImage(*pimage);
  delete pimage;

  ScalarBuffer final_result;
  // only valid on rank 0
  apcomp::ScalarImage final_image = compositor.Composite();
  if(dray::mpi_rank() == 0)
  {
    final_result = convert(final_image, m_field_names);
  }
  return final_result;
#else
  // we have composited locally so there is nothing to do
  return scalar_buffer;
#endif
}

} // namespace dray
