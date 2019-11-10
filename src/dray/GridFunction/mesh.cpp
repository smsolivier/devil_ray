#include <dray/GridFunction/mesh.hpp>
#include <dray/GridFunction/mesh_utils.hpp>
#include <dray/GridFunction/device_mesh.hpp>
#include <dray/dray.hpp>
#include <dray/array_utils.hpp>
#include <dray/aabb.hpp>
#include <dray/point_location.hpp>
#include <dray/policies.hpp>
#include <RAJA/RAJA.hpp>

#include <dray/Element/element.hpp>


namespace dray
{

template <class ElemT>
const BVH Mesh<ElemT>::get_bvh() const
{
  return m_bvh;
}

template<class ElemT>
Mesh<ElemT>::Mesh(const GridFunctionData<3u> &dof_data, int32 poly_order)
  : m_dof_data(dof_data),
    m_poly_order(poly_order)
{
  m_bvh = detail::construct_bvh(*this, m_ref_aabbs);
}



//
//HACK to avoid calling eval_inverse() on 2x3 elements.
//
template <uint32 d>
struct LocateHack { };

// 3D: Works.
template <>
struct LocateHack<3u>
{
  template <class ElemT>
  static bool eval_inverse(
      const ElemT &elem,
      stats::IterativeProfile &iter_prof,
      const Vec<typename ElemT::get_precision,3u> &world_coords,
      const AABB<3u> &guess_domain,
      Vec<typename ElemT::get_precision,3u> &ref_coords,
      bool use_init_guess = false)
  {
    return elem.eval_inverse(iter_prof, world_coords, guess_domain, ref_coords, use_init_guess);
  }

  template <class ElemT>
  static bool eval_inverse(
      const ElemT &elem,
      const Vec<typename ElemT::get_precision,3u> &world_coords,
      const AABB<3u> &guess_domain,
      Vec<typename ElemT::get_precision,3u> &ref_coords,
      bool use_init_guess = false)
  {
    return elem.eval_inverse(world_coords, guess_domain, ref_coords, use_init_guess);
  }
};

// 2D: Dummy, does nothing.
template <>
struct LocateHack<2u>
{
  template <class ElemT>
  static bool eval_inverse(
      const ElemT &elem,
      stats::IterativeProfile &iter_prof,
      const Vec<typename ElemT::get_precision,3u> &world_coords,
      const AABB<2u> &guess_domain,
      Vec<typename ElemT::get_precision,2u> &ref_coords,
      bool use_init_guess = false)
  {
    return false;
  }

  template <class ElemT>
  static bool eval_inverse(
      const ElemT &elem,
      const Vec<typename ElemT::get_precision,3u> &world_coords,
      const AABB<2u> &guess_domain,
      Vec<typename ElemT::get_precision,2u> &ref_coords,
      bool use_init_guess = false)
  {
    return false;
  }
};


template<class ElemT>
Mesh<ElemT>::Mesh()
{
  #warning "need default mesh constructor"
}

template<class ElemT>
AABB<3>
Mesh<ElemT>::get_bounds() const
{
  return m_bvh.m_bounds;
}

template<class ElemT>
template <class StatsType>
void Mesh<ElemT>::locate(Array<int32> &active_idx,
                         Array<Vec<Float,3u>> &wpoints,
                         Array<Location> &locations,
                         StatsType &stats) const
{
  //template <int32 _RefDim>
  //using BShapeOp = BernsteinBasis<T,3>;
  //using ShapeOpType = BShapeOp<3>;

  const int32 size = wpoints.size();
  const int32 size_active = active_idx.size();
  // The results will go in rpoints. Make sure there's room.
  assert((locations.size() >= size_active));

  PointLocator locator(m_bvh);
  //constexpr int32 max_candidates = 5;
  constexpr int32 max_candidates = 100;
  //Size size_active * max_candidates.
  PointLocator::Candidates candidates = locator.locate_candidates(wpoints,
                                                                  active_idx,
                                                                  max_candidates);

  const AABB<dim> *ref_aabb_ptr = m_ref_aabbs.get_device_ptr_const();

  // Initialize outputs to well-defined dummy values.
  Vec<Float,dim> three_point_one_four;
  three_point_one_four = 3.14;

  // Assume that elt_ids and ref_pts are sized to same length as wpoints.
  //assert(elt_ids.size() == ref_pts.size());

  const int32  *active_idx_ptr = active_idx.get_device_ptr_const();

  Location *loc_ptr = locations.get_device_ptr();

  const Vec<Float,3> *wpoints_ptr = wpoints.get_device_ptr_const();
  const int32    *cell_id_ptr = candidates.m_candidates.get_device_ptr_const();
  const int32    *aabb_id_ptr = candidates.m_aabb_ids.get_device_ptr_const();

#ifdef DRAY_STATS
  stats::AppStatsAccess device_appstats = stats.get_device_appstats();

  Array<stats::MattStats> mstats;
  mstats.resize(size);
  stats::MattStats *mstats_ptr = mstats.get_device_ptr();
#endif

  DeviceMesh<ElemT> device_mesh(*this);

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, size_active), [=] DRAY_LAMBDA (int32 aii)
  {
#ifdef DRAY_STATS
    stats::MattStats mstat;
    mstat.construct();
#endif
    const int32 ii = active_idx_ptr[aii];
    Location loc = {-1, {-1.f, -1.f, -1.f}};
    const Vec<Float,3> target_pt = wpoints_ptr[ii];

    // - Use aii to index into candidates.
    // - Use ii to index into wpoints, elt_ids, and ref_pts.

    int32 count = 0;
    int32 el_idx = cell_id_ptr[aii*max_candidates + count];
    int32 aabb_idx = aabb_id_ptr[aii*max_candidates + count];
    Vec<Float,dim> el_coords;
    // For accounting/debugging.
    AABB<> cand_overlap = AABB<>::universe();

    bool found_inside = false;
    int32 steps_taken = 0;
    while(!found_inside && count < max_candidates && el_idx != -1)
    {
      steps_taken = 0;
      const bool use_init_guess = true;

      // For accounting/debugging.
      AABB<> bbox;
      device_mesh.get_elem(el_idx).get_bounds(bbox);
      cand_overlap.intersect(bbox);

      AABB<dim> ref_start_box = ref_aabb_ptr[aabb_idx];

#ifdef DRAY_STATS
      stats::IterativeProfile iter_prof;
      iter_prof.construct();
      mstat.m_candidates++;

      found_inside = LocateHack<ElemT::get_dim()>::template eval_inverse<ElemT>(
          device_mesh.get_elem(el_idx),
          iter_prof,
          target_pt,
          ref_start_box,
          el_coords,
          use_init_guess);

      /// found_inside = device_mesh.get_elem(el_idx).eval_inverse(iter_prof,
      ///                                      target_pt,
      ///                                      ref_start_box,
      ///                                      el_coords,
      ///                                      use_init_guess);  // Much easier than before.
      steps_taken = iter_prof.m_num_iter;
      mstat.m_newton_iters += steps_taken;

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_query_stats_ptr[ii].m_total_tests, 1);

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_query_stats_ptr[ii].m_total_test_iterations,
          steps_taken);

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_elem_stats_ptr[el_idx].m_total_tests, 1);

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_elem_stats_ptr[el_idx].m_total_test_iterations,
          steps_taken);
#else
      found_inside = LocateHack<ElemT::get_dim()>::template eval_inverse<ElemT>(
          device_mesh.get_elem(el_idx),
          target_pt,
          ref_start_box,
          el_coords,
          use_init_guess);

      /// found_inside = device_mesh.get_elem(el_idx).eval_inverse(
      ///                                      target_pt,
      ///                                      ref_start_box,
      ///                                      el_coords,
      ///                                      use_init_guess);
#endif

      if (!found_inside && count < max_candidates-1)
      {
        // Continue searching with the next candidate.
        count++;
        el_idx = cell_id_ptr[aii*max_candidates + count];
        aabb_idx = aabb_id_ptr[aii*max_candidates + count];
      }
    }

    // After testing each candidate, now record the result.
    if (found_inside)
    {
      loc.m_cell_id = el_idx;
      loc.m_ref_pt[0] = el_coords[0];
      loc.m_ref_pt[1] = el_coords[1];
      if(dim == 3)
      {
        loc.m_ref_pt[2] = el_coords[2];
      }

    }

    loc_ptr[ii] = loc;

#ifdef DRAY_STATS

    if (found_inside)
    {
      mstat.m_found = 1;

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_query_stats_ptr[ii].m_total_hits,
          1);

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_query_stats_ptr[ii].m_total_hit_iterations,
          steps_taken);

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_elem_stats_ptr[el_idx].m_total_hits,
          1);

      RAJA::atomicAdd<atomic_policy>(
          &device_appstats.m_elem_stats_ptr[el_idx].m_total_hit_iterations,
          steps_taken);
    }
    mstats_ptr[aii] = mstat;
#endif
  });

#ifdef DRAY_STATS
  stats::StatStore::add_point_stats(wpoints, mstats);
#endif
}


// Explicit instantiations.
//template class MeshAccess<MeshElem<2u, ElemType::Quad, Order::General>>;
//template class MeshAccess<MeshElem<2u, ElemType::Tri, Order::General>>;
//
//template class MeshAccess<MeshElem<3u, ElemType::Quad, Order::General>>;
//template class MeshAccess<MeshElem<3u, ElemType::Tri, Order::General>>;

// Explicit instantiations.
template class Mesh<MeshElem<2u, ElemType::Quad, Order::General>>;
/// template class Mesh<float32, MeshElem<float32, 2u, ElemType::Tri, Order::General>>;
/// template class Mesh<float64, MeshElem<float64, 2u, ElemType::Tri, Order::General>>;

template class Mesh<MeshElem<3u, ElemType::Quad, Order::General>>;
/// template class Mesh<float32, MeshElem<float32, 3u, ElemType::Tri, Order::General>>;   //TODO change ref boxes to SubRef<etype>
/// template class Mesh<float64, MeshElem<float64, 3u, ElemType::Tri, Order::General>>;
}
