#include <dray/filters/mesh_lines.hpp>

#include <dray/array_utils.hpp>
#include <dray/ref_point.hpp>
#include <dray/shaders.hpp>
#include <dray/high_order_intersection.hpp>


namespace dray
{

namespace detail
{

//
// intersect_AABB()
//
// Copied verbatim from triangle_mesh.cpp
//
template <typename T>
DRAY_EXEC_ONLY
bool intersect_AABB(const Vec<float32,4> *bvh,
                    const int32 &currentNode,
                    const Vec<T,3> &orig_dir,
                    const Vec<T,3> &inv_dir,
                    const T& closest_dist,
                    bool &hit_left,
                    bool &hit_right,
                    const T &min_dist) //Find hit after this distance
{
  Vec<float32, 4> first4  = const_get_vec4f(&bvh[currentNode + 0]);
  Vec<float32, 4> second4 = const_get_vec4f(&bvh[currentNode + 1]);
  Vec<float32, 4> third4  = const_get_vec4f(&bvh[currentNode + 2]);
  T xmin0 = first4[0] * inv_dir[0] - orig_dir[0];
  T ymin0 = first4[1] * inv_dir[1] - orig_dir[1];
  T zmin0 = first4[2] * inv_dir[2] - orig_dir[2];
  T xmax0 = first4[3] * inv_dir[0] - orig_dir[0];
  T ymax0 = second4[0] * inv_dir[1] - orig_dir[1];
  T zmax0 = second4[1] * inv_dir[2] - orig_dir[2];
  T min0 = fmaxf(
    fmaxf(fmaxf(fminf(ymin0, ymax0), fminf(xmin0, xmax0)), fminf(zmin0, zmax0)),
    min_dist);
  T max0 = fminf(
    fminf(fminf(fmaxf(ymin0, ymax0), fmaxf(xmin0, xmax0)), fmaxf(zmin0, zmax0)),
    closest_dist);
  hit_left = (max0 >= min0);

  T xmin1 = second4[2] * inv_dir[0] - orig_dir[0];
  T ymin1 = second4[3] * inv_dir[1] - orig_dir[1];
  T zmin1 = third4[0] * inv_dir[2] - orig_dir[2];
  T xmax1 = third4[1] * inv_dir[0] - orig_dir[0];
  T ymax1 = third4[2] * inv_dir[1] - orig_dir[1];
  T zmax1 = third4[3] * inv_dir[2] - orig_dir[2];

  T min1 = fmaxf(
    fmaxf(fmaxf(fminf(ymin1, ymax1), fminf(xmin1, xmax1)), fminf(zmin1, zmax1)),
    min_dist);
  T max1 = fminf(
    fminf(fminf(fmaxf(ymin1, ymax1), fmaxf(xmin1, xmax1)), fmaxf(zmin1, zmax1)),
    closest_dist);
  hit_right = (max1 >= min1);

  return (min0 > min1);
}

//
// candidate_ray_intersection()
//
//   TODO find appropriate place for this function. It is mostly copied from TriangleMesh
//
template <typename T, int32 max_candidates>
Array<int32> candidate_ray_intersection(Array<Ray<T>> rays, const BVH bvh)
{
  const int32 size = rays.size();

  Array<int32> candidates;
  candidates.resize(size * max_candidates);
  array_memset(candidates, -1);

  //const int32 *active_ray_ptr = rays.m_active_rays.get_device_ptr_const();
  const Ray<T> *ray_ptr = rays.get_device_ptr_const();

  const int32 *leaf_ptr = bvh.m_leaf_nodes.get_device_ptr_const();
  const Vec<float32, 4> *inner_ptr = bvh.m_inner_nodes.get_device_ptr_const();

  int32 *candidates_ptr = candidates.get_device_ptr();

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, size), [=] DRAY_LAMBDA (int32 i)
  {

    const Ray<T> &ray = ray_ptr[i];

    T closest_dist = ray.m_far;
    T min_dist = ray.m_near;
    ///int32 hit_idx = -1;
    const Vec<T,3> dir = ray.m_dir;
    Vec<T,3> inv_dir;
    inv_dir[0] = rcp_safe(dir[0]);
    inv_dir[1] = rcp_safe(dir[1]);
    inv_dir[2] = rcp_safe(dir[2]);

    int32 current_node;
    int32 todo[max_candidates];
    int32 stackptr = 0;
    current_node = 0;

    constexpr int32 barrier = -2000000000;
    todo[stackptr] = barrier;

    const Vec<T,3> orig = ray.m_orig;

    Vec<T,3> orig_dir;
    orig_dir[0] = orig[0] * inv_dir[0];
    orig_dir[1] = orig[1] * inv_dir[1];
    orig_dir[2] = orig[2] * inv_dir[2];

    int32 candidate_idx = 0;

    while (current_node != barrier && candidate_idx < max_candidates)
    {
      if (current_node > -1)
      {
        bool hit_left, hit_right;
        bool right_closer = intersect_AABB(inner_ptr,
                                           current_node,
                                           orig_dir,
                                           inv_dir,
                                           closest_dist,
                                           hit_left,
                                           hit_right,
                                           min_dist);

        if (!hit_left && !hit_right)
        {
          current_node = todo[stackptr];
          stackptr--;
        }
        else
        {
          Vec<float32, 4> children = const_get_vec4f(&inner_ptr[current_node + 3]);
          int32 l_child;
          constexpr int32 isize = sizeof(int32);
          memcpy(&l_child, &children[0], isize);
          int32 r_child;
          memcpy(&r_child, &children[1], isize);
          current_node = (hit_left) ? l_child : r_child;

          if (hit_left && hit_right)
          {
            if (right_closer)
            {
              current_node = r_child;
              stackptr++;
              todo[stackptr] = l_child;
            }
            else
            {
              stackptr++;
              todo[stackptr] = r_child;
            }
          }
        }
      } // if inner node

      if (current_node < 0 && current_node != barrier) //check register usage
      {
        current_node = -current_node - 1; //swap the neg address

        // Any leaf bbox we enter is a candidate.
        candidates_ptr[candidate_idx + i * max_candidates] = leaf_ptr[current_node];
        candidate_idx++;

        current_node = todo[stackptr];
        stackptr--;
      } // if leaf node

    } //while

  });

  return candidates;
}

  struct HasCandidate
  {
    int32 m_max_candidates;
    const int32 *m_candidates_ptr;
    DRAY_EXEC bool operator() (int32 ii) const { return (m_candidates_ptr[ii * m_max_candidates] > -1); }
  };

}  // namespace detail

//template Array<RefPoint<float32,3>>
//intersect_mesh_faces(Array<Ray<float32>> rays, const Mesh<float32> &mesh);
//template Array<RefPoint<float64,3>>
//intersect_mesh_faces(Array<Ray<float64>> rays, const Mesh<float64> &mesh);
//
//template Array<Vec<float32,4>>
//mesh_lines<float32>(Array<Ray<float32>> rays, const Mesh<float32> &mesh);
//template Array<Vec<float32,4>>
//mesh_lines<float64>(Array<Ray<float64>> rays, const Mesh<float64> &mesh);

template <typename T>
Array<RefPoint<T,3>> intersect_mesh_faces(Array<Ray<T>> rays, const Mesh<T> &mesh)
{
  constexpr int32 ref_dim = 3;

  // Initialize rpoints to same size as rays, each rpoint set to invalid_refpt.
  Array<RefPoint<T,3>> rpoints;
  rpoints.resize(rays.size());
  const RefPoint<T,ref_dim> invalid_refpt{ -1, {-1,-1,-1} };
  array_memset(rpoints, invalid_refpt);

  // Duplicated from MeshField::intersect_mesh_boundary().

  const Vec<T,3> element_guess = {0.5, 0.5, 0.5};
  const T ray_guess = 1.0;

  // Get intersection candidates for all active rays.
  constexpr int32 max_candidates = 64;
  Array<int32> candidates =
    detail::candidate_ray_intersection<T, max_candidates> (rays, mesh.get_bvh());
  const int32 *candidates_ptr = candidates.get_device_ptr_const();

  const int32 size = rays.size();

    // Define pointers for RAJA kernel.
  MeshAccess<T> device_mesh = mesh.access_device_mesh();
  Ray<T> *ray_ptr = rays.get_device_ptr();
  RefPoint<T,ref_dim> *rpoints_ptr = rpoints.get_device_ptr();

  // For each active ray, loop through candidates until found an intersection.
  RAJA::forall<for_policy>(RAJA::RangeSegment(0, size), [=] DRAY_LAMBDA (const int32 i)
  {
    Ray<T> &ray = ray_ptr[i];
    RefPoint<T,ref_dim> &rpt = rpoints_ptr[i];

    Vec<T,3> ref_coords = element_guess;
    T ray_dist = ray_guess;

    // In case no intersection is found.
    ray.m_active = 0;
    // TODO change comparisons for valid rays to check both near and far.
    ray.m_dist = infinity<T>();

    bool found_inside = false;
    int32 candidate_idx = 0;
    int32 el_idx = candidates_ptr[i*max_candidates + candidate_idx];

    while (!found_inside && candidate_idx < max_candidates && el_idx != -1)
    {
      ref_coords = element_guess;
      ray_dist = ray_guess;
      const bool use_init_guess = true;

#ifdef DRAY_STATS
      stats::IterativeProfile iter_prof;    iter_prof.construct();

      MeshElem<T> mesh_elem = device_mesh.get_elem(el_idx);
      found_inside = Intersector_RayFace<T>::intersect(iter_prof, mesh_elem, ray,
          ref_coords, ray_dist, use_init_guess);
#else
        //TODO after add this to the Intersector_RayFace interface.
        /// found_inside =
        //Intersector_RayFace<T>::intersect(device_mesh.get_elem(el_idx), ray,
        ///     ref_coords, ray_dist, use_init_guess);
#endif
        if (found_inside && ray_dist < ray.m_dist && ray_dist >= ray.m_near)
        {
          rpt.m_el_id = el_idx;
          rpt.m_el_coords = ref_coords;
          ray.m_dist = ray_dist;
        }

        // Continue searching with the next candidate.
        candidate_idx++;
        el_idx = candidates_ptr[i*max_candidates + candidate_idx];

      } // end while

  });  // end RAJA

  return rpoints;
}

template <typename T>
Array<Vec<float32,4>> mesh_lines(Array<Ray<T>> rays, const Mesh<T,3> &mesh)
{
  using Color = Vec<float32,4>;
  constexpr int32 ref_dim = 3;

  Array<Color> color_buffer;
  color_buffer.resize(rays.size());

  // Initialize the color buffer to (0,0,0,0).
  const Color init_color = make_vec4f(0.f, 0.f, 0.f, 0.f);
  const Color bg_color   = make_vec4f(1.f, 1.f, 1.f, 1.f);
  array_memset_vec(color_buffer, init_color);

  // Initialize fragment shader.
  ShadeMeshLines shader;
  const Color face_color = make_vec4f(0.f, 0.f, 0.f, 1.f);
  const Color line_color = make_vec4f(1.f, 1.f, 1.f, 1.f);
  const float32 line_ratio = 0.05;
  shader.set_uniforms(line_color, face_color, line_ratio);

  // Start the rays out at the min distance from calc ray start.
  // Note: Rays that have missed the mesh bounds will have near >= far,
  //       so after the copy, we can detect misses as dist >= far.
  dray::AABB<3> mesh_bounds = mesh.get_bounds();
  calc_ray_start(rays, mesh_bounds);
  Array<int32> active_rays = active_indices(rays);

  // Remove the rays which totally miss the mesh.
  rays = gather(rays, active_rays);

  Array<RefPoint<T,ref_dim>> rpoints = intersect_mesh_faces(rays, mesh);
  Color *color_buffer_ptr = color_buffer.get_device_ptr();
  const RefPoint<T,3> *rpoints_ptr = rpoints.get_device_ptr_const();
  const Ray<T> *rays_ptr = rays.get_device_ptr_const();

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, rpoints.size()), [=] DRAY_LAMBDA (int32 ii)
  {
    const RefPoint<T> &rpt = rpoints_ptr[ii];
    if (rpt.m_el_id != -1)
    {
      Color pixel_color = shader(rpt.m_el_coords);
      color_buffer_ptr[rays_ptr[ii].m_pixel_id] = pixel_color;
    }
  });

  Shader::composite_bg(color_buffer, bg_color);

  return color_buffer;

}

template<typename T>
Array<Vec<float32,4>>
MeshLines::execute(Array<Ray<T>> &rays, DataSet<T> &data_set)
{
  return mesh_lines(rays, data_set.get_mesh());
}

template
Array<Vec<float32,4>>
MeshLines::execute<float64>(Array<Ray<float64>> &rays,
                            DataSet<float64> &data_set);

template
Array<Vec<float32,4>>
MeshLines::execute<float32>(Array<Ray<float32>> &rays,
                            DataSet<float32> &data_set);


};//naemespace dray