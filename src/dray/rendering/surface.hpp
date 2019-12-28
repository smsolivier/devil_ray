// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef DRAY_SURFACE_HPP
#define DRAY_SURFACE_HPP

#include<dray/rendering/traceable.hpp>

namespace dray
{

class Surface : public Traceable
{
protected:
  bool m_draw_mesh;
  float32 m_line_thickness;
  float32 m_sub_res;       // sub resolution of grid lines
public:
  Surface() = delete;
  Surface(DataSet &dataset);
  virtual ~Surface();

  virtual Array<RayHit> nearest_hit(Array<Ray> &rays) override;

  virtual void shade(const Array<Ray> &rays,
                     const Array<RayHit> &hits,
                     const Array<Fragment> &fragments,
                     const Array<PointLight> &lights,
                     Framebuffer &framebuffer) override;

  template<typename MeshElem>
  Array<RayHit> execute(Mesh<MeshElem> &mesh, Array<Ray> &rays);
  void draw_mesh(bool on);
  void line_thickness(const float32 thickness);
};

};//namespace dray

#endif //DRAY_SURFACE_HPP