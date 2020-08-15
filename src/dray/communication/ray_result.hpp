#ifndef DRAY_RAY_RESULT_HPP
#define DRAY_RAY_RESULT_HPP

#include <dray/vec.hpp>

namespace dray
{

struct RayResult
{
  Vec<float32,4> m_color;
  int32 m_pixel_id;
};

} //namespace dray
#endif
