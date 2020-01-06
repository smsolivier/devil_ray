#include <dray/filters/path_lengths.hpp>
#include <dray/uniform_topology.hpp>
#include <dray/policies.hpp>
#include <dray/utils/point_writer.hpp>

namespace dray
{

namespace detail
{

Array<Vec<Float,3>> cell_centers(UniformTopology &topo)
{

  const Vec<int32,3> cell_dims = topo.cell_dims();
  const Vec<Float,3> origin = topo.origin();
  const Vec<Float,3> spacing = topo.spacing();

  const int32 num_cells = cell_dims[0] * cell_dims[1] * cell_dims[2];

  Array<Vec<Float,3>> locations;
  locations.resize(num_cells);
  Vec<Float,3> *loc_ptr = locations.get_device_ptr();

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, num_cells), [=] DRAY_LAMBDA (int32 index)
  {
    Vec<int32,3> cell_id;
    cell_id[0] = index % cell_dims[0];
    cell_id[1] = (index / cell_dims[0]) % cell_dims[1];
    cell_id[2] = index / (cell_dims[0] * cell_dims[1]);

    Vec<Float,3> loc;
    for(int32 i = 0; i < 3; ++i)
    {
      loc[i] = origin[i] + Float(cell_id[i]) * spacing[i] + spacing[i] * 0.5f;
    }

    loc_ptr[index] = loc;
  });

  return locations;
}

struct TraversalState
{
  Vec<Float,3> m_delta_max;
  Vec<Float,3> m_delta;
  Vec<int32,3> m_voxel;
  Vec<Float,3> m_dir;

  // distance to voxel exit from initial point
  DRAY_EXEC
  Float exit() const
  {
    return min(m_delta_max[0], min(m_delta_max[1], m_delta_max[2]));
  }

  // advances to the next voxel along the ray
  DRAY_EXEC void advance()
  {
    int32 advance_dir = 0;
    for(int32 i = 1; i < 3; ++i)
    {
      if(m_delta_max[i] < m_delta_max[advance_dir])
      {
        advance_dir = i;
      }
    }
    m_delta_max[advance_dir] += m_delta[advance_dir];
    m_voxel[advance_dir] += m_dir[advance_dir] < 0.f ? -1 : 1;

    //std::cout<<"Voxel "<<voxel<<"\n";
    //assert(m_voxel[0] >= 0);
    //assert(m_voxel[1] >= 0);
    //assert(m_voxel[2] >= 0);

    //assert(m_voxel[0] < m_dims[0]);
    //assert(m_voxel[1] < m_dims[1]);
    //assert(m_voxel[2] < m_dims[2]);
  }

};

struct DDATraversal
{
  const Vec<int32,3> m_dims;
  const Vec<Float,3> m_origin;
  const Vec<Float,3> m_spacing;

  DDATraversal(UniformTopology &topo)
    : m_dims(topo.cell_dims()),
      m_origin(topo.origin()),
      m_spacing(topo.spacing())
  {

  }

  DRAY_EXEC
  bool is_inside(const Vec<int32, 3>& index) const
  {
    bool inside = true;
    const int32 minIndex = min(index[0], min(index[1], index[2]));
    if(minIndex < 0) inside = false;
    if(index[0] >= m_dims[0]) inside = false;
    if(index[1] >= m_dims[1]) inside = false;
    if(index[2] >= m_dims[2]) inside = false;
    return inside;
  }

  DRAY_EXEC Float
  init_traversal(const Vec<Float,3> &point,
                 const Vec<Float,3> &dir,
                 TraversalState &state) const
  {
    //assert(is_inside(point));
    Vec<Float, 3> temp = point;
    temp = temp - m_origin;
    state.m_voxel[0] = temp[0] / m_spacing[0];
    state.m_voxel[1] = temp[1] / m_spacing[1];
    state.m_voxel[2] = temp[2] / m_spacing[2];
    state.m_dir = dir;

    Vec<Float,3> step;
    step[0] = (dir[0] >= 0.f) ? 1.f : -1.f;
    step[1] = (dir[1] >= 0.f) ? 1.f : -1.f;
    step[2] = (dir[2] >= 0.f) ? 1.f : -1.f;

    Vec<Float,3> next_boundary;
    next_boundary[0] = (Float(state.m_voxel[0]) + step[0]) * m_spacing[0];
    next_boundary[1] = (Float(state.m_voxel[1]) + step[1]) * m_spacing[1];
    next_boundary[2] = (Float(state.m_voxel[2]) + step[2]) * m_spacing[2];

    // correct next boundary for negative directions
    if(step[0] == -1.f) next_boundary[0] += m_spacing[0];
    if(step[1] == -1.f) next_boundary[1] += m_spacing[1];
    if(step[2] == -1.f) next_boundary[2] += m_spacing[2];

    // distance to next voxel boundary
    state.m_delta_max[0] = (dir[0] != 0.f) ?
      (next_boundary[0] - (point[0] - m_origin[0])) / dir[0] : infinity<Float>();

    state.m_delta_max[1] = (dir[1] != 0.f) ?
      (next_boundary[1] - (point[1] - m_origin[1])) / dir[1] : infinity<Float>();

    state.m_delta_max[2] = (dir[2] != 0.f) ?
      (next_boundary[2] - (point[2] - m_origin[2])) / dir[2] : infinity<Float>();

    // distance along ray to traverse x,y, and z of a voxel
    state.m_delta[0] = (dir[0] != 0) ? m_spacing[0] / dir[0] * step[0] : infinity<Float>();
    state.m_delta[1] = (dir[1] != 0) ? m_spacing[1] / dir[1] * step[1] : infinity<Float>();
    state.m_delta[2] = (dir[2] != 0) ? m_spacing[2] / dir[2] * step[2] : infinity<Float>();

    Vec<Float,3> exit_boundary;
    exit_boundary[0] = step[0] < 0.f ? 0.f : Float(m_dims[0]) * m_spacing[0];
    exit_boundary[1] = step[1] < 0.f ? 0.f : Float(m_dims[1]) * m_spacing[1];
    exit_boundary[2] = step[2] < 0.f ? 0.f : Float(m_dims[2]) * m_spacing[2];

    if(step[0] == -1.f) exit_boundary[0] += m_spacing[0];
    if(step[1] == -1.f) exit_boundary[1] += m_spacing[1];
    if(step[2] == -1.f) exit_boundary[2] += m_spacing[2];

    Vec<Float,3> exit_dist;
    // distance to grid exit
    exit_dist[0] = (dir[0] != 0.f) ?
      (exit_boundary[0] - (point[0] - m_origin[0])) / dir[0] : infinity<Float>();

    exit_dist[1] = (dir[1] != 0.f) ?
      (exit_boundary[1] - (point[1] - m_origin[1])) / dir[1] : infinity<Float>();

    exit_dist[2] = (dir[2] != 0.f) ?
      (exit_boundary[2] - (point[2] - m_origin[2])) / dir[2] : infinity<Float>();

    //std::cout<<"Init voxel "<<voxel<<"\n";

    return min(exit_dist[0], min(exit_dist[1], exit_dist[2]));
  }



};


} // namespace detail

PathLengths::PathLengths()
 :  m_x_res(10),
    m_y_res(10),
    m_width(2.f),
    m_height(2.f),
    m_origin({{0.f, 0.f, 15.f}}),
    m_normal({{0.f, 0.f, 1.f}}),
    m_x_dir({{1.f, 0.f, 0.f}})
{
}

Array<Vec<Float,3>>
PathLengths::generate_pixels()
{
  const Float pixel_width = Float(m_x_res) / m_width;
  const Float pixel_height = Float(m_y_res) / m_height;

  // better be orthogonal
  assert(dot(m_normal, m_x_dir) == 0.f);

  Vec<Float,3> y_dir = cross(m_normal, m_x_dir);

  // avoid lambda capturing 'this'
  Vec<Float,3> x_dir = m_x_dir;
  const int32 width = m_x_res;

  Vec<Float,3> start;
  start = m_origin;
  start += (pixel_width / 2.f) * x_dir;
  start += (pixel_height / 2.f) * y_dir;

  const int32 num_pixels = m_x_res * m_y_res;

  Array<Vec<Float,3>> pixels;
  pixels.resize(num_pixels);
  Vec<Float,3> *pixel_ptr = pixels.get_device_ptr();

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, num_pixels), [=] DRAY_LAMBDA (int32 index)
  {
    const int32 x = index % width;
    const int32 y = index / width;

    Vec<Float,3> loc;
    loc = start + (x * pixel_width * x_dir) + (y * pixel_height * y_dir);
    pixel_ptr[index] = loc;
  });

  return pixels;
}

void go(Array<Vec<Float,3>> &pixels, UniformTopology &topo)
{
  const detail::DDATraversal dda(topo);
  const Vec<Float,3> *pixel_ptr = pixels.get_device_ptr_const();
  const int32 size = pixels.size();

  RAJA::forall<for_policy>(RAJA::RangeSegment(0, size), [=] DRAY_LAMBDA (int32 index)
  {
    Vec<Float,3> pixel = pixel_ptr[index];
    Vec<Float,3> loc = {{1.5f, 0.5f, 0.5f}}; // dummy point
    Vec<Float,3> dir = pixel - loc;
    dir.normalize();
    detail::TraversalState state;
    dda.init_traversal(loc, dir, state);

    Float distance = 0.f;

    while(dda.is_inside(state.m_voxel))
    {
      const Float voxel_exit = state.exit();
      const Float length = voxel_exit - distance;
      if(index == 1) std::cout<<state.m_voxel<<" length "<<length<<"\n";
      // do stuff
      distance = voxel_exit;
      state.advance();
    }

  });
}

void PathLengths::execute(DataSet &data_set)
{
  Array<Vec<Float,3>> pixels = generate_pixels();
  write_points(pixels);

  TopologyBase *topo = data_set.topology();
  if(dynamic_cast<UniformTopology*>(topo) != nullptr)
  {
    std::cout<<"Boom\n";
    UniformTopology *uni_topo = dynamic_cast<UniformTopology*>(topo);
    go(pixels, *uni_topo);

  }
}

};//namespace dray

