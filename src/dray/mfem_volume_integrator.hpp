#ifndef DRAY_MFEM_VOLUME_INTEGRATOR_HPP
#define DRAY_MFEM_VOLUME_INTEGRATOR_HPP

#include <dray/mfem_mesh.hpp>
#include <dray/mfem_grid_function.hpp>

namespace dray
{

class MFEMVolumeIntegrator
{
protected:
  MFEMMeshField  m_mesh;
  float32        m_sample_dist;

  MFEMVolumeIntegrator(); 
public:
  MFEMVolumeIntegrator(MFEMMeshField &mesh); 
  ~MFEMVolumeIntegrator(); 
  
  template<typename T>
  void            integrate(Ray<T> &rays);
  
};

} // namespace dray

#endif
