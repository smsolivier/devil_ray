#include <dray/rendering/volume_sampler.hpp>

#include <dray/dray.hpp>

#ifdef DRAY_MPI_ENABLED
#include <mpi.h>
#endif

namespace dray
{

VolumeSampler::VolumeSampler(Collection &collection)
  : m_collection(collection)
{
  init();
}

VolumeSampler::~VolumeSampler()
{
}

void
VolumeSampler::input(Collection &collection)
{
  m_collection = collection;
  init();
}

void
VolumeSampler::field(const std::string field)
{
  m_field = field;
  m_field_range = m_collection.range(m_field);
}

std::string
VolumeSampler::field() const
{
  return m_field;
}

void VolumeSampler::init()
{
  m_domains.clear();
  m_bounds_map.clear();

  int32 num_domains = m_collection.size();

  int domain_offset = 0;
#ifdef DRAY_MPI_ENABLED
  int comm_size = dray::mpi_size();
  int rank = dray::mpi_rank();
  MPI_Comm comm = MPI_Comm_f2c(dray::mpi_comm());

  int *domains_per_rank = new int[comm_size];
  MPI_Allgather(&num_domains, 1, MPI_INT, domains_per_rank, 1, MPI_INT, comm);
  for(int i = 0; i < rank; ++i)
  {
    domain_offset += domains_per_rank[i];
  }
  delete[] domains_per_rank;
#endif

  for(int i = 0; i < num_domains; ++i)
  {
    DataSet domain = m_collection.domain(i);
    AABB<3> aabb = domain.topology()->bounds();
    m_bounds_map.add_block(i + domain_offset, aabb);

    DomainSampler sampler(domain);
    m_domains[i+domain_offset] = sampler;
  }

  m_bounds_map.build();


}


} // namespace dray
