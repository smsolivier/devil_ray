#include <dray/bounds_map.hpp>
#include <dray/error.hpp>
#include <dray/dray.hpp>
#include <dray/linear_bvh_builder.hpp>

#ifdef DRAY_MPI_ENABLED
#include <mpi.h>
#endif

namespace dray
{

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

void BoundsMap::build()
{
  int size = m_bounds.size();
#if DRAY_MPI_ENABLEDV
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

  int idx = 0;
  for(auto bounds : m_bounds)
  {
    dom_ids_ptr[idx] = bounds.first;
    aabbs_ptr[idx] = bounds.second;
    idx++;
  }

  LinearBVHBuilder builder;
  m_bvh = builder.construct(aabbs, domain_ids);
}


} // namespace dray
