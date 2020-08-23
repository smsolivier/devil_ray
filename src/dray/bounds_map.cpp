#include <dray/bounds_map.hpp>
#include <dray/error.hpp>
#include <dray/dray.hpp>
#include <dray/linear_bvh_builder.hpp>


#ifdef DRAY_MPI_ENABLED
#include <mpi.h>
#endif

namespace dray
{

namespace detail
{

bool intersect_AABB_dist_host(const Vec<float32,4> *bvh,
                              const int32 &currentNode,
                              const Vec<Float,3> &orig_dir,
                              const Vec<Float,3> &inv_dir,
                              const Float& closest_dist,
                              bool &hit_left,
                              bool &hit_right,
                              Float &ldist,
                              Float &rdist,
                              const Float &min_dist) //Find hit after this distance
{
  Vec<float32, 4> first4  = bvh[currentNode + 0];
  Vec<float32, 4> second4 = bvh[currentNode + 1];
  Vec<float32, 4> third4  = bvh[currentNode + 2];
  Float xmin0 = first4[0] * inv_dir[0] - orig_dir[0];
  Float ymin0 = first4[1] * inv_dir[1] - orig_dir[1];
  Float zmin0 = first4[2] * inv_dir[2] - orig_dir[2];
  Float xmax0 = first4[3] * inv_dir[0] - orig_dir[0];
  Float ymax0 = second4[0] * inv_dir[1] - orig_dir[1];
  Float zmax0 = second4[1] * inv_dir[2] - orig_dir[2];
  Float min0 = fmaxf(
    fmaxf(fmaxf(fminf(ymin0, ymax0), fminf(xmin0, xmax0)), fminf(zmin0, zmax0)),
    min_dist);
  Float max0 = fminf(
    fminf(fminf(fmaxf(ymin0, ymax0), fmaxf(xmin0, xmax0)), fmaxf(zmin0, zmax0)),
    closest_dist);
  hit_left = (max0 >= min0);

  Float xmin1 = second4[2] * inv_dir[0] - orig_dir[0];
  Float ymin1 = second4[3] * inv_dir[1] - orig_dir[1];
  Float zmin1 = third4[0] * inv_dir[2] - orig_dir[2];
  Float xmax1 = third4[1] * inv_dir[0] - orig_dir[0];
  Float ymax1 = third4[2] * inv_dir[1] - orig_dir[1];
  Float zmax1 = third4[3] * inv_dir[2] - orig_dir[2];

  Float min1 = fmaxf(
    fmaxf(fmaxf(fminf(ymin1, ymax1), fminf(xmin1, xmax1)), fminf(zmin1, zmax1)),
    min_dist);
  Float max1 = fminf(
    fminf(fminf(fmaxf(ymin1, ymax1), fmaxf(xmin1, xmax1)), fmaxf(zmin1, zmax1)),
    closest_dist);
  hit_right = (max1 >= min1);

  ldist = min0;
  rdist = min1;

  return (min0 > min1);
}

void intersect_ray_host(Ray &ray, BVH &bvh, std::vector<int32> &hits)
{

  const int32 *leaf_ptr = bvh.m_leaf_nodes.get_host_ptr_const();
  const int32 *aabb_ids_ptr = bvh.m_aabb_ids.get_host_ptr_const();
  const Vec<float32, 4> *inner_ptr = bvh.m_inner_nodes.get_host_ptr_const();

  Float closest_dist = ray.m_far;
  Float min_dist = ray.m_near;
  const Vec<Float,3> dir = ray.m_dir;
  Vec<Float,3> inv_dir;
  inv_dir[0] = rcp_safe(dir[0]);
  inv_dir[1] = rcp_safe(dir[1]);
  inv_dir[2] = rcp_safe(dir[2]);

  int32 current_node;
  int32 todo[64];
  int32 stackptr = 0;
  current_node = 0;

  constexpr int32 barrier = -2000000000;
  todo[stackptr] = barrier;

  const Vec<Float,3> orig = ray.m_orig;

  Vec<Float,3> orig_dir;
  orig_dir[0] = orig[0] * inv_dir[0];
  orig_dir[1] = orig[1] * inv_dir[1];
  orig_dir[2] = orig[2] * inv_dir[2];

  while (current_node != barrier)
  {
    if (current_node > -1)
    {
      Float ldist, rdist;
      bool hit_left, hit_right;
      bool right_closer = intersect_AABB_dist_host(inner_ptr,
                                                   current_node,
                                                   orig_dir,
                                                   inv_dir,
                                                   closest_dist,
                                                   hit_left,
                                                   hit_right,
                                                   ldist,
                                                   rdist,
                                                   min_dist);
      if (!hit_left && !hit_right)
      {
        current_node = todo[stackptr];
        stackptr--;
      }
      else
      {
        Vec<float32, 4> children = inner_ptr[current_node + 3];
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

    if (current_node < 0 && current_node != barrier)
    {
      //if(current_distance > closest_dist) std::cout<<"B";
      current_node = -current_node - 1; //swap the neg address
      hits.push_back(current_node);

      current_node = todo[stackptr];
      stackptr--;
    } // if leaf node
  } //while
}

} // namespace detail

void BoundsMap::clear()
{
  m_bounds.clear();
  m_rank_map.clear();
}

void BoundsMap::add_block(int32 id, const AABB<3> &bounds)
{
  if (m_bounds.find(id) == m_bounds.end())
  {
    m_bounds[id] = bounds;
  }
  else
  {
    DRAY_ERROR("Duplicate block id "<<id);
  }

  m_rank_map[id] = dray::mpi_rank();
}

int BoundsMap::get_rank(const int32 &block_id)
{
  auto it = m_rank_map.find(block_id);
  int32 rank = -1;
  if(it != m_rank_map.end())
  {
    rank = m_rank_map[block_id];
  }
  return rank;
}

void BoundsMap::find(Ray &ray, std::vector<int32> &domain_ids)
{
  domain_ids.clear();
  detail::intersect_ray_host(ray, m_bvh, domain_ids);
}

void BoundsMap::build()
{
  int size = m_bounds.size();
#if DRAY_MPI_ENABLED
  int rank;
  int procs;

  MPI_Comm comm = MPI_Comm_f2c(dray::mpi_comm());
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &procs);

  int *dom_counts = new int[procs];
  int *box_counts = new int[procs];

  MPI_Allgather(&size, 1, MPI_INT, dom_counts, 1, MPI_INT, comm);

  if(rank == 0)
  {
    for(int i = 0; i < procs; ++i)
    {
      std::cout<<"rank "<<i<<" has "<<dom_counts[i]<<" boxes\n";
    }
  }

  // prefix sum to build incoming buffers offsets
  int *box_offsets = new int[procs];
  int *dom_offsets = new int[procs];
  box_offsets[0] = 0;
  dom_offsets[0] = 0;
  box_counts[0] = dom_counts[0] * 6;
  for(int i = 1; i < procs; ++i)
  {
    box_offsets[i] = box_offsets[i-1] + dom_counts[i-1] * 6;
    dom_offsets[i] = dom_offsets[i-1] + dom_counts[i-1];
    box_counts[i] = dom_counts[i] * 6;
  }

  int total_boxs = dom_offsets[procs - 1] + dom_counts[procs - 1];
  double *box_send_buff = new double[size * 6];
  int *dom_send_buff = new int[size];

  int counter = 0;
  for (auto it = m_bounds.begin(); it != m_bounds.end(); it++)
  {
     const int offset = counter * 6;
     box_send_buff[offset + 0] = it->second.min()[0];
     box_send_buff[offset + 1] = it->second.max()[0];
     box_send_buff[offset + 2] = it->second.min()[1];
     box_send_buff[offset + 3] = it->second.max()[1];
     box_send_buff[offset + 4] = it->second.min()[2];
     box_send_buff[offset + 5] = it->second.max()[2];
     dom_send_buff[counter] = it->first;
     counter++;
  }
  //std::cout<<"total boxs "<<total_boxs<<"\n";
  double *box_rec_buff = new double[total_boxs * 6];
  int *dom_rec_buff = new int[total_boxs];

  MPI_Allgatherv(box_send_buff,
                 size*6,
                 MPI_DOUBLE,
                 box_rec_buff,
                 box_counts,
                 box_offsets,
                 MPI_DOUBLE,
                 comm);

  MPI_Allgatherv(dom_send_buff,
                 size,
                 MPI_INT,
                 dom_rec_buff,
                 dom_counts,
                 dom_offsets,
                 MPI_INT,
                 comm);

  // repopulate with global information
  clear();

  //build a map of rank that handles empty counts
  int *rank_map = new int[total_boxs];
  int idx = 0;
  for(int i = 0; i < procs; ++i)
  {
    for(int d = 0; d < dom_counts[i]; ++d)
    {
      rank_map[idx] = i;
      ++idx;
    }
  }

  for(int i = 0; i < total_boxs; ++i)
  {
    const int offset = i * 6;
    int dom_id = dom_rec_buff[i];
    AABB<3> &bounds = m_bounds[dom_id];
    Vec<Float,3> mins,maxs;
    mins[0] = box_rec_buff[offset + 0];
    maxs[0] = box_rec_buff[offset + 1];
    mins[1] = box_rec_buff[offset + 2];
    maxs[1] = box_rec_buff[offset + 3];
    mins[2] = box_rec_buff[offset + 4];
    maxs[2] = box_rec_buff[offset + 5];
    bounds.include(mins);
    bounds.include(maxs);

    m_rank_map[dom_id] = rank_map[i];
    if(rank == 0)
    {
      std::cout<<"domain "<<dom_id<<" rank "<<rank_map[i]<<"\n";
    }
  }

  delete[] dom_send_buff;
  delete[] dom_rec_buff;
  delete[] box_send_buff;
  delete[] box_rec_buff;
  delete[] dom_offsets;
  delete[] box_offsets;
  delete[] dom_counts;
  delete[] box_counts;
  delete[] rank_map;
#endif

  Array<AABB<3>> aabbs;
  aabbs.resize(m_bounds.size());

  Array<int32> domain_ids;
  domain_ids.resize(m_bounds.size());

  AABB<3> *aabbs_ptr = aabbs.get_host_ptr();
  int32 *dom_ids_ptr = domain_ids.get_host_ptr();

  int itr = 0;
  for(auto bounds : m_bounds)
  {
    dom_ids_ptr[itr] = bounds.first;
    aabbs_ptr[itr] = bounds.second;
    itr++;
  }

  LinearBVHBuilder builder;
  m_bvh = builder.construct(aabbs, domain_ids);
}


} // namespace dray
