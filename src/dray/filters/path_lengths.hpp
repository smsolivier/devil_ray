#ifndef DRAY_PATH_LENGTHS_HPP
#define DRAY_PATH_LENGTHS_HPP

#include <dray/data_set.hpp>
#include <dray/Element/elem_utils.hpp>

namespace dray
{

class PathLengths
{
protected:
  int32 m_x_res; // detector x resolution
  int32 m_y_res; // detector y resolution
  Float m_width; // detector width
  Float m_height;// detector height
  Vec<Float,3> m_origin; // position of bottom quad (bottom left)
  Vec<Float,3> m_normal; // quad orientation
  Vec<Float,3> m_x_dir;  // quad roll about the normal
public:
  PathLengths();
  void execute(DataSet &data_set);
  Array<Vec<Float,3>> generate_pixels();
  //template<class ElemT>
  //DataSet execute(Mesh<ElemT> &mesh, DataSet &data_set);
};

};//namespace dray

#endif//DRAY_PATH_LENGTHS_HPP
