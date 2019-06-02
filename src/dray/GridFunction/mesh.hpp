#ifndef DRAY_MESH_HPP
#define DRAY_MESH_HPP

#include <dray/GridFunction/grid_function_data.hpp>
#include <dray/Element/element.hpp>
#include <dray/newton_solver.hpp>
#include <dray/vec.hpp>
#include <dray/exports.hpp>

namespace dray
{

  /*
   * @class MeshElem
   * @brief Hexahedral element with (p+1)^3 dofs representing a transformation in the Bernstein basis.
   */
  template <typename T, int32 dim = 3>
  class MeshElem : public Element<T,dim,dim>
  {
    public:
    using Element<T,dim,dim>::construct;
    DRAY_EXEC static MeshElem create(int32 el_id, int32 poly_order, const int32 *ctrl_idx_ptr, const Vec<T,dim> *val_ptr);

    /* Forward evaluation: See Element::eval()   (which for now is ElTransOp::eval(). */

    //
    // eval_inverse() : Try to locate the point in reference space. Return false if not contained.
    DRAY_EXEC bool eval_inverse(const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess = false) const;
    DRAY_EXEC bool eval_inverse(stats::IterativeProfile &iter_prof, const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess = false) const;
  };


  /*
   * @class MeshAccess
   * @brief Device-safe access to a collection of elements (just knows about the geometry, not fields).
   */
  template <typename T, int32 dim = 3>
  struct MeshAccess
  {
    const int32 *m_idx_ptr;
    const Vec<T,dim> *m_val_ptr;
    const int32  m_poly_order;

    //
    // get_elem()
    DRAY_EXEC MeshElem<T,dim> get_elem(int32 el_idx) const;

    // world2ref()
    DRAY_EXEC bool world2ref(int32 el_idx, const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess = false) const;
    DRAY_EXEC bool world2ref(stats::IterativeProfile &iter_prof, int32 el_idx, const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess = false) const;
  };
  

  /*
   * @class Mesh
   * @brief Host-side access to a collection of elements (just knows about the geometry, not fields).
   */
  template <typename T, int32 dim = 3>
  class Mesh
  {
    public:
      Mesh() = delete;  // For now, probably need later.
      Mesh(const GridFunctionData<T,dim> &dof_data, int32 poly_order) : m_dof_data(dof_data), m_poly_order(poly_order) {}
      
      //
      // access_device_mesh() : Must call this BEFORE capture to RAJA lambda.
      MeshAccess<T,dim> access_device_mesh() const;

      //
      // access_host_mesh()
      MeshAccess<T,dim> access_host_mesh() const;

      //
      // get_poly_order()
      int32 get_poly_order() { return m_poly_order; }

      //
      // get_num_elem()
      int32 get_num_elem() { return m_dof_data.get_num_elem(); }

      //
      // get_dof_data()  // TODO should this be removed?
      GridFunctionData<T,dim> get_dof_data() { return m_dof_data; }

    protected:
      GridFunctionData<T,dim> m_dof_data;
      int32 m_poly_order;
  };

}







// Implementations (could go in a .tcc file and include that at the bottom of .hpp)

namespace dray
{

  // ---------------- //
  // MeshElem methods //
  // ---------------- //

  template <typename T, int32 dim>
  DRAY_EXEC MeshElem<T,dim> MeshElem<T,dim>::create(int32 el_id, int32 poly_order, const int32 *ctrl_idx_ptr, const Vec<T,dim> *val_ptr)
  {
    MeshElem<T,dim> ret;
    ret.construct(el_id, poly_order, ctrl_idx_ptr, val_ptr);
    return ret;
  }


  template <typename T, int32 dim>
  DRAY_EXEC bool MeshElem<T,dim>::eval_inverse(const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess) const
  {
    stats::IterativeProfile iter_prof;   iter_prof.construct();
    return eval_inverse(iter_prof, world_coords, ref_coords, use_init_guess);
  }

  template <typename T, int32 dim>
  DRAY_EXEC bool MeshElem<T,dim>::eval_inverse(stats::IterativeProfile &iter_prof, const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess) const
  {
    if (!use_init_guess)
      for (int32 d = 0; d < dim; d++)
        ref_coords[d] = 0.5;

    using IterativeMethod = IterativeMethod<T>;

    // Newton step to solve inverse of geometric transformation.
    struct Stepper {
      DRAY_EXEC typename IterativeMethod::StepStatus operator()
        (Vec<T,dim> &x) const
      {
        Vec<T,dim> delta_y;
        Vec<Vec<T,dim>,dim> j_col;
        Matrix<T,dim,dim> jacobian;
        m_transf.eval(x, delta_y, j_col);
        delta_y = m_target - delta_y;

        for (int32 rdim = 0; rdim < dim; rdim++)
          jacobian.set_col(rdim, j_col[rdim]);

        bool inverse_valid;
        Vec<T,dim> delta_x;
        delta_x = matrix_mult_inv(jacobian, delta_y, inverse_valid);

        if (!inverse_valid)
          return IterativeMethod::Abort;

        x = x + delta_x;
        return IterativeMethod::Continue;
      }

      Element<T,dim,dim> m_transf;
      Vec<T,dim> m_target;

    } stepper{ *this, world_coords};

    //TODO somewhere else in the program, figure out how to set the precision
    //based on the gradient and the image resolution.
    const T tol_ref = 1e-6;
    const int32 max_steps = 100;

    // Find solution.
    return (IterativeMethod::solve(iter_prof, stepper, ref_coords, max_steps, tol_ref) == IterativeMethod::Converged
      && Element<T,dim,dim>::is_inside(ref_coords));
  }


  // ------------------ //
  // MeshAccess methods //
  // ------------------ //

  //
  // get_elem()
  template <typename T, int32 dim>
  DRAY_EXEC MeshElem<T,dim> MeshAccess<T,dim>::get_elem(int32 el_idx) const
  {
    // We are just going to assume that the elements in the data store
    // are in the same position as their id, el_id==el_idx.
    MeshElem<T,dim> ret;
    ret.construct(el_idx, m_poly_order, m_idx_ptr, m_val_ptr);
    return ret;
  }

  //
  // world2ref()
  template <typename T, int32 dim>
  DRAY_EXEC bool MeshAccess<T,dim>::world2ref(int32 el_idx, const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess) const
  {
    return get_elem(el_idx).eval_inverse(world_coords, ref_coords, use_init_guess);
  }
  template <typename T, int32 dim>
  DRAY_EXEC bool MeshAccess<T,dim>::world2ref(stats::IterativeProfile &iter_prof, int32 el_idx, const Vec<T,dim> &world_coords, Vec<T,dim> &ref_coords, bool use_init_guess) const
  {
    return get_elem(el_idx).eval_inverse(iter_prof, world_coords, ref_coords, use_init_guess);
  }


  // ---------------- //
  // Mesh methods     //
  // ---------------- //

  //
  // access_device_mesh()
  template <typename T, int32 dim>
  MeshAccess<T,dim> Mesh<T,dim>::access_device_mesh() const
  {
    return { m_dof_data.m_ctrl_idx.get_device_ptr_const(),
             m_dof_data.m_values.get_device_ptr_const(), m_poly_order };
  }

  //
  // access_host_mesh()
  template <typename T, int32 dim>
  MeshAccess<T,dim> Mesh<T,dim>::access_host_mesh() const
  {
    return { m_dof_data.m_ctrl_idx.get_host_ptr_const(),
             m_dof_data.m_values.get_host_ptr_const(),
             m_poly_order };
  }

} // namespace dray


#endif//DRAY_MESH_HPP