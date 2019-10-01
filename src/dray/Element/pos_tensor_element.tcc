#ifndef DRAY_POS_TENSOR_ELEMENT_HPP
#define DRAY_POS_TENSOR_ELEMENT_HPP

/**
 * @file pos_tensor_element.hpp
 * @brief Partial template specialization of Element_impl
 *        for tensor (i.e. hex and quad) elements.
 */

#include <dray/Element/element.hpp>
#include <dray/integer_utils.hpp>   // MultinomialCoeff
#include <dray/vec.hpp>

#include <dray/Element/bernstein_basis.hpp>  // get_sub_coefficient

namespace dray
{
  //
  // QuadElement_impl
  //
  template <typename T, uint32 dim, uint32 ncomp, int32 P>
  using QuadElement_impl = Element_impl<T, dim, ncomp, ElemType::Quad, P>;


  template <typename T, uint32 dim>
  class QuadRefSpace
  {
    public:
      DRAY_EXEC static bool is_inside(const Vec<T,dim> &ref_coords);  //TODO
      DRAY_EXEC static void clamp_to_domain(Vec<T,dim> &ref_coords);  //TODO
      DRAY_EXEC static Vec<T,dim> project_to_domain(const Vec<T,dim> &r1, const Vec<T,dim> &r2);  //TODO
  };


  // Specialize SubRef for Quad type.
  template <>
  struct ElemTypeAttributes<ElemType::Quad>
  {
    template <uint32 dim>
    using SubRef = AABB<dim>;
  };

  // TODO add get_sub_bounds to each specialization.

  // ---------------------------------------------------------------------------


  // -----
  // Implementations
  // -----

  template <typename T, uint32 dim>
  DRAY_EXEC bool QuadRefSpace<T,dim>::is_inside(const Vec<T,dim> &ref_coords)
  {
    T min_val = 2.f;
    T max_val = -1.f;
    for (int32 d = 0; d < dim; d++)
    {
      min_val = min(ref_coords[d], min_val);
      max_val = max(ref_coords[d], max_val);
    }
    return (min_val >= 0.f - epsilon<T>()) && (max_val <= 1.f + epsilon<T>());
  }
  
  template <typename T, uint32 dim>
  DRAY_EXEC void QuadRefSpace<T,dim>::clamp_to_domain(Vec<T,dim> &ref_coords)
  {
    //TODO
  }
  
  template <typename T, uint32 dim>
  DRAY_EXEC Vec<T,dim> QuadRefSpace<T,dim>::project_to_domain(const Vec<T,dim> &r1, const Vec<T,dim> &r2)
  {
    return {0.0}; //TODO
  }



  // ---------------------------------------------------------------------------

  // Template specialization (Quad type, general order).
  //
  // Assume dim <= 3.
  //
  template <typename T, uint32 dim, uint32 ncomp>
  class Element_impl<T, dim, ncomp, ElemType::Quad, Order::General> : public QuadRefSpace<T,dim>
  {
    protected:
      SharedDofPtr<Vec<T, ncomp>> m_dof_ptr;
      uint32 m_order;
    public:
      DRAY_EXEC void construct(SharedDofPtr<Vec<T, ncomp>> dof_ptr, int32 poly_order)
      {
        m_dof_ptr = dof_ptr;
        m_order = poly_order;
      }
      DRAY_EXEC int32 get_order() const    { return m_order; }
      DRAY_EXEC int32 get_num_dofs() const { return get_num_dofs(m_order); }
      DRAY_EXEC static constexpr int32 get_num_dofs(int32 order) { return intPow(order+1, dim); }

      DRAY_EXEC Vec<T, ncomp> eval(const Vec<T,dim> &r) const
      {
        using DofT = Vec<T, ncomp>;
        using PtrT = SharedDofPtr<Vec<T, ncomp>>;

        //TODO
        DofT answer; answer = 0;
        return answer;
      }

      DRAY_EXEC Vec<T, ncomp> eval_d( const Vec<T,dim> &ref_coords,
                                      Vec<Vec<T,ncomp>,dim> &out_derivs) const
      {
        // Directly evaluate a Bernstein polynomial with a hybrid of Horner's rule and accumulation of powers:
        //     V = 0.0;  xpow = 1.0;
        //     for(i)
        //     {
        //       V = V*(1-x) + C[i]*xpow*nchoosek(p,i);
        //       xpow *= x;
        //     }
        //
        // Indirectly evaluate a high-order Bernstein polynomial, by directly evaluating
        // the two parent lower-order Bernstein polynomials, and mixing with weights {(1-x), x}.
        //
        // Indirectly evaluate the derivative of a high-order Bernstein polynomial, by directly
        // evaluating the two parent lower-order Bernstein polynomials, and mixing with weights {-p, p}.

        using DofT = Vec<T, ncomp>;
        using PtrT = SharedDofPtr<Vec<T, ncomp>>;

        DofT zero;
        zero = 0;

        const T u = (dim > 0 ? ref_coords[0] : 0.0);
        const T v = (dim > 1 ? ref_coords[1] : 0.0);
        const T w = (dim > 2 ? ref_coords[2] : 0.0);
        const T ubar = 1.0 - u;
        const T vbar = 1.0 - v;
        const T wbar = 1.0 - w;

        const int32 p1 = (dim >= 1 ? m_order : 0);
        const int32 p2 = (dim >= 2 ? m_order : 0);
        const int32 p3 = (dim >= 3 ? m_order : 0);

        int32 B[MaxPolyOrder];
        if (m_order >= 1)
        {
          BinomialCoeff binomial_coeff;
          binomial_coeff.construct(m_order - 1);
          for (int32 ii = 0; ii <= m_order - 1; ii++)
          {
            B[ii] = binomial_coeff.get_val();
            binomial_coeff.slide_over(0);
          }
        }

        int32 cidx = 0;  // Index into element dof indexing space.

        // Compute and combine order (p-1) values to get order (p) values/derivatives.
        // https://en.wikipedia.org/wiki/Bernstein_polynomial#Properties

        DofT val_u, val_v, val_w;
        DofT        deriv_u;
        Vec<DofT,2> deriv_uv;
        Vec<DofT,3> deriv_uvw;

        // Level3 set up.
        T wpow = 1.0;
        Vec<DofT,3> val_w_L, val_w_R;  // Second/third columns are derivatives in lower level.
        val_w_L = zero;
        val_w_R = zero;
        for (int32 ii = 0; ii <= p3; ii++)
        {
          // Level2 set up.
          T vpow = 1.0;
          Vec<DofT,2> val_v_L, val_v_R;  // Second column is derivative in lower level.
          val_v_L = zero;
          val_v_R = zero;

          for (int32 jj = 0; jj <= p2; jj++)
          {
            // Level1 set up.
            T upow = 1.0;
            DofT val_u_L = zero, val_u_R = zero;           // L and R can be combined --> val, deriv.
            DofT C = m_dof_ptr[cidx++];
            for (int32 kk = 1; kk <= p1; kk++)
            {
              // Level1 accumulation.
              val_u_L = val_u_L * ubar + C * (B[kk-1] * upow);
              C = m_dof_ptr[cidx++];
              val_u_R = val_u_R * ubar + C * (B[kk-1] * upow);
              upow *= u;
            }//kk

            // Level1 result.
            val_u = (p1 > 0 ? val_u_L * ubar + val_u_R * u : C);
            deriv_u = (val_u_R - val_u_L) * p1;

            // Level2 accumulation.
            if (jj > 0)
            {
              val_v_R[0] = val_v_R[0] * vbar + val_u   * (B[jj-1] * vpow);
              val_v_R[1] = val_v_R[1] * vbar + deriv_u * (B[jj-1] * vpow);
              vpow *= v;
            }
            if (jj < p2)
            {
              val_v_L[0] = val_v_L[0] * vbar + val_u   * (B[jj] * vpow);
              val_v_L[1] = val_v_L[1] * vbar + deriv_u * (B[jj] * vpow);
            }
          }//jj

          // Level2 result.
          val_v       = (p2 > 0 ? val_v_L[0] * vbar + val_v_R[0] * v : val_u);
          deriv_uv[0] = (p2 > 0 ? val_v_L[1] * vbar + val_v_R[1] * v : deriv_u);
          deriv_uv[1] = (val_v_R[0] - val_v_L[0]) * p2;

          // Level3 accumulation.
          if (ii > 0)
          {
            val_w_R[0] = val_w_R[0] * wbar + val_v       * (B[ii-1] * wpow);
            val_w_R[1] = val_w_R[1] * wbar + deriv_uv[0] * (B[ii-1] * wpow);
            val_w_R[2] = val_w_R[2] * wbar + deriv_uv[1] * (B[ii-1] * wpow);
            wpow *= w;
          }
          if (ii < p3)
          {
            val_w_L[0] = val_w_L[0] * wbar + val_v       * (B[ii] * wpow);
            val_w_L[1] = val_w_L[1] * wbar + deriv_uv[0] * (B[ii] * wpow);
            val_w_L[2] = val_w_L[2] * wbar + deriv_uv[1] * (B[ii] * wpow);
          }
        }//ii

        // Level3 result.
        val_w        = (p3 > 0 ? val_w_L[0] * wbar + val_w_R[0] * w : val_v);
        deriv_uvw[0] = (p3 > 0 ? val_w_L[1] * wbar + val_w_R[1] * w : deriv_uv[0]);
        deriv_uvw[1] = (p3 > 0 ? val_w_L[2] * wbar + val_w_R[2] * w : deriv_uv[1]);
        deriv_uvw[2] = (val_w_R[0] - val_w_L[0]) * p3;

        if (dim > 0) out_derivs[0] = deriv_uvw[0];
        if (dim > 1) out_derivs[1] = deriv_uvw[1];
        if (dim > 2) out_derivs[2] = deriv_uvw[2];

        return val_w;
      }

      DRAY_EXEC void get_sub_bounds(const AABB<dim> &sub_ref, AABB<ncomp> &aabb) const;
  };


  //
  // get_sub_bounds()
  template <typename T, uint32 dim, uint32 ncomp>
  DRAY_EXEC void
  Element_impl<T, dim, ncomp, ElemType::Quad, Order::General>::get_sub_bounds(
      const AABB<dim> &sub_ref,
      AABB<ncomp> &aabb) const
  {
    // Initialize.
    aabb.reset();

    const int32 num_dofs = get_num_dofs();

    using DofT = Vec<T, ncomp>;
    using PtrT = SharedDofPtr<Vec<T, ncomp>>;

    if (m_order <= 3) // TODO find the optimal threshold, if there is one.
    {
      // Get the sub-coefficients all at once in a block.
      switch(m_order)
      {
        case 1:
        {
          constexpr int32 POrder = 1;
          MultiVec<T, dim, ncomp, POrder> sub_nodes = sub_element_fixed_order<T, dim, ncomp, POrder, PtrT>(sub_ref.m_ranges, m_dof_ptr);
          for (int32 ii = 0; ii < num_dofs; ii++)
            aabb.include(sub_nodes.linear_idx(ii));
        }
        break;

        case 2:
        {
          constexpr int32 POrder = 2;
          MultiVec<T, dim, ncomp, POrder> sub_nodes = sub_element_fixed_order<T, dim, ncomp, POrder, PtrT>(sub_ref.m_ranges, m_dof_ptr);
          for (int32 ii = 0; ii < num_dofs; ii++)
            aabb.include(sub_nodes.linear_idx(ii));
        }
        break;

        case 3:
        {
          constexpr int32 POrder = 3;
          MultiVec<T, dim, ncomp, POrder> sub_nodes = sub_element_fixed_order<T, dim, ncomp, POrder, PtrT>(sub_ref.m_ranges, m_dof_ptr);
          for (int32 ii = 0; ii < num_dofs; ii++)
            aabb.include(sub_nodes.linear_idx(ii));
        }
        break;
      }
    }
    else
    {
      // Get each sub-coefficient one at a time.
      for (int32 i0 = 0; i0 <= (dim >= 1 ? m_order : 0); i0++)
        for (int32 i1 = 0; i1 <= (dim >= 2 ? m_order : 0); i1++)
          for (int32 i2 = 0; i2 <= (dim >= 3 ? m_order : 0); i2++)
          {
            Vec<T,ncomp> sub_node =
              // TODO move out of bernstein_basis.hpp
              BernsteinBasis<T,dim>::template get_sub_coefficient<PtrT, ncomp>(sub_ref.m_ranges,
                                                                               m_dof_ptr,
                                                                               m_order,
                                                                               i0,
                                                                               i1,
                                                                               i2);
              aabb.include(sub_node);
          }
    }
  }






  // ---------------------------------------------------------------------------


  // Template specialization (Tensor type, 0th order).
  //
  template <typename T, uint32 dim, uint32 ncomp>
  class Element_impl<T, dim, ncomp, ElemType::Quad, Order::Constant> : public QuadRefSpace<T,dim>
  {
    protected:
      SharedDofPtr<Vec<T, ncomp>> m_dof_ptr;
    public:
      DRAY_EXEC void construct(SharedDofPtr<Vec<T, ncomp>> dof_ptr, int32 p) { m_dof_ptr = dof_ptr; }
      DRAY_EXEC static constexpr int32 get_order() { return 0; }
      DRAY_EXEC static constexpr int32 get_num_dofs() { return 1; }
      DRAY_EXEC static constexpr int32 get_num_dofs(int32) { return 1; }

      // Get value without derivative.
      DRAY_EXEC Vec<T, ncomp> eval(const Vec<T,dim> &ref_coords) const
      {
        return *m_dof_ptr;
      }

      // Get value with derivative.
      DRAY_EXEC Vec<T, ncomp> eval_d( const Vec<T,dim> &ref_coords,
                                      Vec<Vec<T,ncomp>,dim> &out_derivs) const
      {
        for (int d = 0; d < dim; d++)
          out_derivs[d] = 0;

        return *m_dof_ptr;
      }
  };


  // Template specialization (Quad type, 1st order, 2D).
  //
  template <typename T, uint32 ncomp>
  class Element_impl<T, 2u, ncomp, ElemType::Quad, Order::Linear> : public QuadRefSpace<T,2u>
  {
    protected:
      SharedDofPtr<Vec<T, ncomp>> m_dof_ptr;
    public:
      DRAY_EXEC void construct(SharedDofPtr<Vec<T, ncomp>> dof_ptr, int32 p) { m_dof_ptr = dof_ptr; }
      DRAY_EXEC static constexpr int32 get_order() { return 1; }
      DRAY_EXEC static constexpr int32 get_num_dofs() { return 4; }
      DRAY_EXEC static constexpr int32 get_num_dofs(int32) { return 4; }

      // Get value without derivative.
      DRAY_EXEC Vec<T, ncomp> eval(const Vec<T,2u> &r) const
      {
        return m_dof_ptr[0] * (1-r[0]) * (1-r[1]) +
               m_dof_ptr[1] *    r[0]  * (1-r[1]) +
               m_dof_ptr[2] * (1-r[0]) *    r[1]  +
               m_dof_ptr[3] *    r[0]  *    r[1];
      }

      // Get value with derivative.
      DRAY_EXEC Vec<T, ncomp> eval_d( const Vec<T,2u> &r,
                                      Vec<Vec<T,ncomp>,2u> &out_derivs) const
      {
        out_derivs[0] = (m_dof_ptr[1] - m_dof_ptr[0]) * (1-r[1]) +
                        (m_dof_ptr[3] - m_dof_ptr[2]) *    r[1];

        out_derivs[1] = (m_dof_ptr[2] - m_dof_ptr[0]) * (1-r[0]) +
                        (m_dof_ptr[3] - m_dof_ptr[1]) *    r[0];

        return m_dof_ptr[0] * (1-r[0]) * (1-r[1]) +
               m_dof_ptr[1] *    r[0]  * (1-r[1]) +
               m_dof_ptr[2] * (1-r[0]) *    r[1]  +
               m_dof_ptr[3] *    r[0]  *    r[1];
      }
  };


  // Template specialization (Quad type, 1st order, 3D).
  //
  template <typename T, uint32 ncomp>
  class Element_impl<T, 3u, ncomp, ElemType::Quad, Order::Linear> : public QuadRefSpace<T,3u>
  {
    protected:
      SharedDofPtr<Vec<T, ncomp>> m_dof_ptr;
    public:
      DRAY_EXEC void construct(SharedDofPtr<Vec<T, ncomp>> dof_ptr, int32 p) { m_dof_ptr = dof_ptr; }
      DRAY_EXEC static constexpr int32 get_order() { return 1; }
      DRAY_EXEC static constexpr int32 get_num_dofs() { return 8; }
      DRAY_EXEC static constexpr int32 get_num_dofs(int32) { return 8; }

      // Get value without derivative.
      DRAY_EXEC Vec<T, ncomp> eval(const Vec<T,3u> &r) const
      {
        return m_dof_ptr[0] * (1-r[0]) * (1-r[1]) * (1-r[2]) +
               m_dof_ptr[1] *    r[0]  * (1-r[1]) * (1-r[2]) +
               m_dof_ptr[2] * (1-r[0]) *    r[1]  * (1-r[2]) +
               m_dof_ptr[3] *    r[0]  *    r[1]  * (1-r[2]) +
               m_dof_ptr[4] * (1-r[0]) * (1-r[1]) *    r[2]  +
               m_dof_ptr[5] *    r[0]  * (1-r[1]) *    r[2]  +
               m_dof_ptr[6] * (1-r[0]) *    r[1]  *    r[2]  +
               m_dof_ptr[7] *    r[0]  *    r[1]  *    r[2]  ;
      }

      // Get value with derivative.
      DRAY_EXEC Vec<T, ncomp> eval_d( const Vec<T,3u> &r,
                                    Vec<Vec<T,ncomp>,3u> &out_derivs) const
      {
        out_derivs[0] = (m_dof_ptr[1] - m_dof_ptr[0]) * (1-r[1]) * (1-r[2]) +
                        (m_dof_ptr[3] - m_dof_ptr[2]) *    r[1]  * (1-r[2]) +
                        (m_dof_ptr[5] - m_dof_ptr[4]) * (1-r[1]) *    r[2]  +
                        (m_dof_ptr[7] - m_dof_ptr[6]) *    r[1]  *    r[2]  ;

        out_derivs[1] = (m_dof_ptr[2] - m_dof_ptr[0]) * (1-r[0]) * (1-r[2]) +
                        (m_dof_ptr[3] - m_dof_ptr[1]) *    r[0]  * (1-r[2]) +
                        (m_dof_ptr[6] - m_dof_ptr[4]) * (1-r[0]) *    r[2]  +
                        (m_dof_ptr[7] - m_dof_ptr[5]) *    r[0]  *    r[2]  ;

        out_derivs[2] = (m_dof_ptr[4] - m_dof_ptr[0]) * (1-r[0]) * (1-r[1]) +
                        (m_dof_ptr[5] - m_dof_ptr[1]) *    r[0]  * (1-r[1]) +
                        (m_dof_ptr[6] - m_dof_ptr[2]) * (1-r[0]) *    r[1]  +
                        (m_dof_ptr[7] - m_dof_ptr[3]) *    r[0]  *    r[1]  ;

        return m_dof_ptr[0] * (1-r[0]) * (1-r[1]) * (1-r[2]) +
               m_dof_ptr[1] *    r[0]  * (1-r[1]) * (1-r[2]) +
               m_dof_ptr[2] * (1-r[0]) *    r[1]  * (1-r[2]) +
               m_dof_ptr[3] *    r[0]  *    r[1]  * (1-r[2]) +
               m_dof_ptr[4] * (1-r[0]) * (1-r[1]) *    r[2]  +
               m_dof_ptr[5] *    r[0]  * (1-r[1]) *    r[2]  +
               m_dof_ptr[6] * (1-r[0]) *    r[1]  *    r[2]  +
               m_dof_ptr[7] *    r[0]  *    r[1]  *    r[2]  ;
      }
  };





  // Template specialization (Quad type, 2nd order, 2D).
  //
  template <typename T, uint32 ncomp>
  class Element_impl<T, 2u, ncomp, ElemType::Quad, Order::Quadratic> : public QuadRefSpace<T,2u>
  {
    protected:
      SharedDofPtr<Vec<T, ncomp>> m_dof_ptr;
    public:
      DRAY_EXEC void construct(SharedDofPtr<Vec<T, ncomp>> dof_ptr, int32 p) { m_dof_ptr = dof_ptr; }
      DRAY_EXEC static constexpr int32 get_order() { return 2; }
      DRAY_EXEC static constexpr int32 get_num_dofs() { return IntPow<3,2u>::val; }
      DRAY_EXEC static constexpr int32 get_num_dofs(int32) { return IntPow<3,2u>::val; }

      DRAY_EXEC Vec<T, ncomp> eval(const Vec<T,2u> &r) const
      {
        // Shape functions. Quadratic has 3 1D shape functions on each axis.
        T su[3] = { (1-r[0])*(1-r[0]),  2*r[0]*(1-r[0]),  r[0]*r[0] };
        T sv[3] = { (1-r[1])*(1-r[1]),  2*r[1]*(1-r[1]),  r[1]*r[1] };

        return m_dof_ptr[0] * su[0] * sv[0] +
               m_dof_ptr[1] * su[1] * sv[0] +
               m_dof_ptr[2] * su[2] * sv[0] +
               m_dof_ptr[3] * su[0] * sv[1] +
               m_dof_ptr[4] * su[1] * sv[1] +
               m_dof_ptr[5] * su[2] * sv[1] +
               m_dof_ptr[6] * su[0] * sv[2] +
               m_dof_ptr[7] * su[1] * sv[2] +
               m_dof_ptr[8] * su[2] * sv[2] ;
      }

      DRAY_EXEC Vec<T, ncomp> eval_d( const Vec<T,2u> &r,
                                      Vec<Vec<T,ncomp>,2u> &out_derivs) const
      {
        // Shape functions. Quadratic has 3 1D shape functions on each axis.
        T su[3] = { (1-r[0])*(1-r[0]),  2*r[0]*(1-r[0]),  r[0]*r[0] };
        T sv[3] = { (1-r[1])*(1-r[1]),  2*r[1]*(1-r[1]),  r[1]*r[1] };

        // Shape derivatives.
        T dsu[3] = { -1+r[0],  1-r[0]-r[0],  r[0] };
        T dsv[3] = { -1+r[1],  1-r[1]-r[1],  r[1] };

        out_derivs[0] = m_dof_ptr[0] * dsu[0] * sv[0] +
                        m_dof_ptr[1] * dsu[1] * sv[0] +
                        m_dof_ptr[2] * dsu[2] * sv[0] +
                        m_dof_ptr[3] * dsu[0] * sv[1] +
                        m_dof_ptr[4] * dsu[1] * sv[1] +
                        m_dof_ptr[5] * dsu[2] * sv[1] +
                        m_dof_ptr[6] * dsu[0] * sv[2] +
                        m_dof_ptr[7] * dsu[1] * sv[2] +
                        m_dof_ptr[8] * dsu[2] * sv[2] ;

        out_derivs[1] = m_dof_ptr[0] * su[0] * dsv[0] +
                        m_dof_ptr[1] * su[1] * dsv[0] +
                        m_dof_ptr[2] * su[2] * dsv[0] +
                        m_dof_ptr[3] * su[0] * dsv[1] +
                        m_dof_ptr[4] * su[1] * dsv[1] +
                        m_dof_ptr[5] * su[2] * dsv[1] +
                        m_dof_ptr[6] * su[0] * dsv[2] +
                        m_dof_ptr[7] * su[1] * dsv[2] +
                        m_dof_ptr[8] * su[2] * dsv[2] ;

        return m_dof_ptr[0] * su[0] * sv[0] +
               m_dof_ptr[1] * su[1] * sv[0] +
               m_dof_ptr[2] * su[2] * sv[0] +
               m_dof_ptr[3] * su[0] * sv[1] +
               m_dof_ptr[4] * su[1] * sv[1] +
               m_dof_ptr[5] * su[2] * sv[1] +
               m_dof_ptr[6] * su[0] * sv[2] +
               m_dof_ptr[7] * su[1] * sv[2] +
               m_dof_ptr[8] * su[2] * sv[2] ;
      }
  };


  // Template specialization (Quad type, 2nd order, 3D).
  //
  template <typename T, uint32 ncomp>
  class Element_impl<T, 3u, ncomp, ElemType::Quad, Order::Quadratic> : public QuadRefSpace<T,3u>
  {
    protected:
      SharedDofPtr<Vec<T, ncomp>> m_dof_ptr;
    public:
      DRAY_EXEC void construct(SharedDofPtr<Vec<T, ncomp>> dof_ptr, int32 p) { m_dof_ptr = dof_ptr; }
      DRAY_EXEC static constexpr int32 get_order() { return 2; }
      DRAY_EXEC static constexpr int32 get_num_dofs() { return IntPow<3,3u>::val; }
      DRAY_EXEC static constexpr int32 get_num_dofs(int32) { return IntPow<3,3u>::val; }

      DRAY_EXEC Vec<T, ncomp> eval(const Vec<T,3u> &r) const
      {
        //TODO
        Vec<T, ncomp> answer; answer = 0;
        return answer;
      }

      DRAY_EXEC Vec<T, ncomp> eval_d( const Vec<T,3u> &r,
                                      Vec<Vec<T,ncomp>,3u> &out_derivs) const
      {
        // Shape functions. Quadratic has 3 1D shape functions on each axis.
        T su[3] = { (1-r[0])*(1-r[0]),  2*r[0]*(1-r[0]),  r[0]*r[0] };
        T sv[3] = { (1-r[1])*(1-r[1]),  2*r[1]*(1-r[1]),  r[1]*r[1] };
        T sw[3] = { (1-r[2])*(1-r[2]),  2*r[2]*(1-r[2]),  r[2]*r[2] };

        // Shape derivatives.
        T dsu[3] = { -1+r[0],  1-r[0]-r[0],  r[0] };
        T dsv[3] = { -1+r[1],  1-r[1]-r[1],  r[1] };
        T dsw[3] = { -1+r[2],  1-r[2]-r[2],  r[2] };

        out_derivs[0] =  m_dof_ptr[0]  * dsu[0] * sv[0] * sw[0] +
                         m_dof_ptr[1]  * dsu[1] * sv[0] * sw[0] +
                         m_dof_ptr[2]  * dsu[2] * sv[0] * sw[0] +
                         m_dof_ptr[3]  * dsu[0] * sv[1] * sw[0] +
                         m_dof_ptr[4]  * dsu[1] * sv[1] * sw[0] +
                         m_dof_ptr[5]  * dsu[2] * sv[1] * sw[0] +
                         m_dof_ptr[6]  * dsu[0] * sv[2] * sw[0] +
                         m_dof_ptr[7]  * dsu[1] * sv[2] * sw[0] +
                         m_dof_ptr[8]  * dsu[2] * sv[2] * sw[0] +

                         m_dof_ptr[9]  * dsu[0] * sv[0] * sw[1] +
                         m_dof_ptr[10] * dsu[1] * sv[0] * sw[1] +
                         m_dof_ptr[11] * dsu[2] * sv[0] * sw[1] +
                         m_dof_ptr[12] * dsu[0] * sv[1] * sw[1] +
                         m_dof_ptr[13] * dsu[1] * sv[1] * sw[1] +
                         m_dof_ptr[14] * dsu[2] * sv[1] * sw[1] +
                         m_dof_ptr[15] * dsu[0] * sv[2] * sw[1] +
                         m_dof_ptr[16] * dsu[1] * sv[2] * sw[1] +
                         m_dof_ptr[17] * dsu[2] * sv[2] * sw[1] +

                         m_dof_ptr[18] * dsu[0] * sv[0] * sw[2] +
                         m_dof_ptr[19] * dsu[1] * sv[0] * sw[2] +
                         m_dof_ptr[20] * dsu[2] * sv[0] * sw[2] +
                         m_dof_ptr[21] * dsu[0] * sv[1] * sw[2] +
                         m_dof_ptr[22] * dsu[1] * sv[1] * sw[2] +
                         m_dof_ptr[23] * dsu[2] * sv[1] * sw[2] +
                         m_dof_ptr[24] * dsu[0] * sv[2] * sw[2] +
                         m_dof_ptr[25] * dsu[1] * sv[2] * sw[2] +
                         m_dof_ptr[26] * dsu[2] * sv[2] * sw[2] ;

        out_derivs[1] =  m_dof_ptr[0]  * su[0] * dsv[0] * sw[0] +
                         m_dof_ptr[1]  * su[1] * dsv[0] * sw[0] +
                         m_dof_ptr[2]  * su[2] * dsv[0] * sw[0] +
                         m_dof_ptr[3]  * su[0] * dsv[1] * sw[0] +
                         m_dof_ptr[4]  * su[1] * dsv[1] * sw[0] +
                         m_dof_ptr[5]  * su[2] * dsv[1] * sw[0] +
                         m_dof_ptr[6]  * su[0] * dsv[2] * sw[0] +
                         m_dof_ptr[7]  * su[1] * dsv[2] * sw[0] +
                         m_dof_ptr[8]  * su[2] * dsv[2] * sw[0] +

                         m_dof_ptr[9]  * su[0] * dsv[0] * sw[1] +
                         m_dof_ptr[10] * su[1] * dsv[0] * sw[1] +
                         m_dof_ptr[11] * su[2] * dsv[0] * sw[1] +
                         m_dof_ptr[12] * su[0] * dsv[1] * sw[1] +
                         m_dof_ptr[13] * su[1] * dsv[1] * sw[1] +
                         m_dof_ptr[14] * su[2] * dsv[1] * sw[1] +
                         m_dof_ptr[15] * su[0] * dsv[2] * sw[1] +
                         m_dof_ptr[16] * su[1] * dsv[2] * sw[1] +
                         m_dof_ptr[17] * su[2] * dsv[2] * sw[1] +

                         m_dof_ptr[18] * su[0] * dsv[0] * sw[2] +
                         m_dof_ptr[19] * su[1] * dsv[0] * sw[2] +
                         m_dof_ptr[20] * su[2] * dsv[0] * sw[2] +
                         m_dof_ptr[21] * su[0] * dsv[1] * sw[2] +
                         m_dof_ptr[22] * su[1] * dsv[1] * sw[2] +
                         m_dof_ptr[23] * su[2] * dsv[1] * sw[2] +
                         m_dof_ptr[24] * su[0] * dsv[2] * sw[2] +
                         m_dof_ptr[25] * su[1] * dsv[2] * sw[2] +
                         m_dof_ptr[26] * su[2] * dsv[2] * sw[2] ;

        out_derivs[2] =  m_dof_ptr[0]  * su[0] * sv[0] * dsw[0] +
                         m_dof_ptr[1]  * su[1] * sv[0] * dsw[0] +
                         m_dof_ptr[2]  * su[2] * sv[0] * dsw[0] +
                         m_dof_ptr[3]  * su[0] * sv[1] * dsw[0] +
                         m_dof_ptr[4]  * su[1] * sv[1] * dsw[0] +
                         m_dof_ptr[5]  * su[2] * sv[1] * dsw[0] +
                         m_dof_ptr[6]  * su[0] * sv[2] * dsw[0] +
                         m_dof_ptr[7]  * su[1] * sv[2] * dsw[0] +
                         m_dof_ptr[8]  * su[2] * sv[2] * dsw[0] +

                         m_dof_ptr[9]  * su[0] * sv[0] * dsw[1] +
                         m_dof_ptr[10] * su[1] * sv[0] * dsw[1] +
                         m_dof_ptr[11] * su[2] * sv[0] * dsw[1] +
                         m_dof_ptr[12] * su[0] * sv[1] * dsw[1] +
                         m_dof_ptr[13] * su[1] * sv[1] * dsw[1] +
                         m_dof_ptr[14] * su[2] * sv[1] * dsw[1] +
                         m_dof_ptr[15] * su[0] * sv[2] * dsw[1] +
                         m_dof_ptr[16] * su[1] * sv[2] * dsw[1] +
                         m_dof_ptr[17] * su[2] * sv[2] * dsw[1] +

                         m_dof_ptr[18] * su[0] * sv[0] * dsw[2] +
                         m_dof_ptr[19] * su[1] * sv[0] * dsw[2] +
                         m_dof_ptr[20] * su[2] * sv[0] * dsw[2] +
                         m_dof_ptr[21] * su[0] * sv[1] * dsw[2] +
                         m_dof_ptr[22] * su[1] * sv[1] * dsw[2] +
                         m_dof_ptr[23] * su[2] * sv[1] * dsw[2] +
                         m_dof_ptr[24] * su[0] * sv[2] * dsw[2] +
                         m_dof_ptr[25] * su[1] * sv[2] * dsw[2] +
                         m_dof_ptr[26] * su[2] * sv[2] * dsw[2] ;

        return m_dof_ptr[0]  * su[0] * sv[0] * sw[0] +
               m_dof_ptr[1]  * su[1] * sv[0] * sw[0] +
               m_dof_ptr[2]  * su[2] * sv[0] * sw[0] +
               m_dof_ptr[3]  * su[0] * sv[1] * sw[0] +
               m_dof_ptr[4]  * su[1] * sv[1] * sw[0] +
               m_dof_ptr[5]  * su[2] * sv[1] * sw[0] +
               m_dof_ptr[6]  * su[0] * sv[2] * sw[0] +
               m_dof_ptr[7]  * su[1] * sv[2] * sw[0] +
               m_dof_ptr[8]  * su[2] * sv[2] * sw[0] +

               m_dof_ptr[9]  * su[0] * sv[0] * sw[1] +
               m_dof_ptr[10] * su[1] * sv[0] * sw[1] +
               m_dof_ptr[11] * su[2] * sv[0] * sw[1] +
               m_dof_ptr[12] * su[0] * sv[1] * sw[1] +
               m_dof_ptr[13] * su[1] * sv[1] * sw[1] +
               m_dof_ptr[14] * su[2] * sv[1] * sw[1] +
               m_dof_ptr[15] * su[0] * sv[2] * sw[1] +
               m_dof_ptr[16] * su[1] * sv[2] * sw[1] +
               m_dof_ptr[17] * su[2] * sv[2] * sw[1] +

               m_dof_ptr[18] * su[0] * sv[0] * sw[2] +
               m_dof_ptr[19] * su[1] * sv[0] * sw[2] +
               m_dof_ptr[20] * su[2] * sv[0] * sw[2] +
               m_dof_ptr[21] * su[0] * sv[1] * sw[2] +
               m_dof_ptr[22] * su[1] * sv[1] * sw[2] +
               m_dof_ptr[23] * su[2] * sv[1] * sw[2] +
               m_dof_ptr[24] * su[0] * sv[2] * sw[2] +
               m_dof_ptr[25] * su[1] * sv[2] * sw[2] +
               m_dof_ptr[26] * su[2] * sv[2] * sw[2] ;
      }
  };



  // Template specialization (Quad type, 3rd order).
  //
  template <typename T, uint32 dim, uint32 ncomp>
  class Element_impl<T, dim, ncomp, ElemType::Quad, Order::Cubic> : public QuadRefSpace<T,dim>
  {
    protected:
      SharedDofPtr<Vec<T, ncomp>> m_dof_ptr;
    public:
      DRAY_EXEC void construct(SharedDofPtr<Vec<T, ncomp>> dof_ptr, int32 p) { m_dof_ptr = dof_ptr; }
      DRAY_EXEC static constexpr int32 get_order() { return 3; }
      DRAY_EXEC static constexpr int32 get_num_dofs() { return IntPow<4,dim>::val; }
      DRAY_EXEC static constexpr int32 get_num_dofs(int32) { return IntPow<4,dim>::val; }

      DRAY_EXEC Vec<T, ncomp> eval(const Vec<T,dim> &r) const
      {
        //TODO
        Vec<T, ncomp> answer; answer = 0;
        return answer;
      }

      DRAY_EXEC Vec<T, ncomp> eval_d( const Vec<T,dim> &ref_coords,
                                      Vec<Vec<T,ncomp>,dim> &out_derivs) const
      {
        //TODO
        Vec<T, ncomp> answer; answer = 0;
        return answer;
      }
  };




}//namespace dray

#endif// DRAY_POS_TENSOR_ELEMENT_HPP
