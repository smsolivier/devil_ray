// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef DRAY_ELEM_OPS_HPP
#define DRAY_ELEM_OPS_HPP

#include <dray/types.hpp>
#include <dray/Element/elem_attr.hpp>
#include <dray/Element/element.hpp>
#include <dray/Element/dof_access.hpp>

namespace dray
{

  namespace detail {
    constexpr int32 cartesian_to_tri_idx(int32 i, int32 j, int32 elen)
    {
      // i runs fastest, j slowest.
      // There are a total of (elen)(elen+1)/2 vertices in the triangle.
      // (idx - i) counts the number of vertices below the cap, so
      //
      //   (elen)(elen+1)/2 - (idx - i) = (elen-j)(elen-j+1)/2
      //
      //   j(1 + 2*elen - j)/2 + i = idx

      return (2*elen + 1 - j)*j/2 + i;
    }

    constexpr int32 cartesian_to_tet_idx(int32 i, int32 j, int32 k, int32 e)
    {
      // i runs fastest, k slowest.
      // There are a total of (elen)(elen+1)(elen+2)/6 vertices in the tetrahedron.
      // (idx - cartesian_to_tri_idx(i,j,elen-k)) counts
      // the number of vertices below the cap, so
      //
      //   (elen)(elen+1)(elen+2)/6 - (idx - (2*elen + 1 - j)*j/2 - i)
      //   = (elen-k)(elen-k+1)(elen-k+2)/6
      //
      //   (e)(e+1)(e+2)/6 - (e-k)(e+1-k)(e+2-k)/6 + (2e + 1 - j)*j/2 + i = idx
      //
      //   ((k - 3e - 3)(k) + (3e + 6)e + 2)k/6 + (2e + 1 - j)*j/2 + i = idx

      return (((-1-e)*3+k)*k + (3*e + 6)*e + 2)*k/6 + (2*e + 1 - j)*j/2 + i;
    }
  }


  /// // get_sub_bounds<Simplex>
  /// template <int32 dim, int32 ncomp, int32 P>
  /// DRAY_EXEC void get_sub_bounds(
  ///     const ShapeTAG,                            // Tag for shape
  ///     const OrderPolicy<P> order_p,              // Tag/data for order policy
  ///     WriteDofPtr<Vec<Float, ncomp>> &dof_ptr,                 // dofs read and written here
  ///     const Split<ElemType::Simplex> &split)
  /// {
  ///   //TODO split the triangle element and use coefficients from split element.
  ///   //For now it just uses the non-split coefficients
  ///   //and returns bounds for entire element.

  ///   const int num_dofs = get_num_dofs(ShapeTAG{}, order_p);

  ///
  /// }

  namespace eops
  {

    template <int32 P>
    struct HexEdgeWalker
    {
      HexEdgeWalker(const OrderPolicy<P> order_p, const int32 eid)
        : m_order_p(order_p),
          m_p(eattr::get_order(order_p)),
          m_edge_base( (m_p+1)*(m_p+1)*(m_p*hex_props::hex_eoffset2(eid))
                             + (m_p+1)*(m_p*hex_props::hex_eoffset1(eid))
                                   + 1*(m_p*hex_props::hex_eoffset0(eid)) ),
          m_edge_stride( hex_props::hex_estride(eid, m_p+1) )
      {}

      const OrderPolicy<P> m_order_p;
      const int32 m_p;
      const int32 m_edge_base;
      const int32 m_edge_stride;
    };

    template <int32 P>
    struct HexFaceWalker
    {
      HexFaceWalker(const OrderPolicy<P> order_p, const int32 fid)
        : m_order_p(order_p),
          m_p(eattr::get_order(order_p)),
          m_f_base( (m_p+1)*(m_p+1)*(m_p*hex_props::hex_foffset2(fid))
                           +(m_p+1)*(m_p*hex_props::hex_foffset1(fid))
                                   +(m_p*hex_props::hex_foffset0(fid)) ),
          m_f_stride0( hex_props::hex_fstrideU(fid, m_p+1) ),
          m_f_stride1( hex_props::hex_fstrideV(fid, m_p+1) )
      {}

      const OrderPolicy<P> m_order_p;
      const int32 m_p;
      const int32 m_f_base;
      const int32 m_f_stride0;
      const int32 m_f_stride1;
    };


    /** eval_d_edge(ShapeHex, Linear) */
    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d_edge( ShapeHex,
                                             const OrderPolicy<Linear> order_p,
                                             const int32 eid,
                                             const ReadDofPtr<Vec<Float, ncomp>> &C,
                                             const Vec<Float, 1> &rc,
                                             Vec<Vec<Float, ncomp>, 1> &out_deriv )
    {
      constexpr int32 p = eattr::get_order(order_p.as_cxp());
      const HexEdgeWalker<Linear> hew(order_p.as_cxp(), eid);
      const Vec<Float, ncomp> C0 = C[ hew.m_edge_base + 0 * hew.m_edge_stride ];
      const Vec<Float, ncomp> C1 = C[ hew.m_edge_base + 1 * hew.m_edge_stride ];
      out_deriv[0] = (C1 - C0)*p;
      return (C1 * rc[0]) + (C0 * (1-rc[0]));
    }

    /** eval_d_edge(ShapeHex, Quadratic) */
    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d_edge( ShapeHex,
                                             const OrderPolicy<Quadratic> order_p,
                                             const int32 eid,
                                             const ReadDofPtr<Vec<Float, ncomp>> &C,
                                             const Vec<Float, 1> &rc,
                                             Vec<Vec<Float, ncomp>, 1> &out_deriv )
    {
      constexpr int32 p = eattr::get_order(order_p.as_cxp());
      const HexEdgeWalker<Quadratic> hew(order_p.as_cxp(), eid);
      Vec<Float, ncomp> C0 = C[ hew.m_edge_base + 0 * hew.m_edge_stride ];
      Vec<Float, ncomp> C1 = C[ hew.m_edge_base + 1 * hew.m_edge_stride ];
      Vec<Float, ncomp> C2 = C[ hew.m_edge_base + 2 * hew.m_edge_stride ];

      C0 = (C1 * rc[0]) + (C0 * (1-rc[0]));
      C1 = (C2 * rc[0]) + (C1 * (1-rc[0]));

      out_deriv[0] = (C1 - C0)*p;
      return (C1 * rc[0]) + (C0 * (1-rc[0]));
    }

    struct BinomialCoeffTable
    {
      // TODO specify gpu 'constant memory' for binomial coefficients.
      int32 m_table[MaxPolyOrder+1];

      BinomialCoeffTable(int32 p)
      {
        BinomialCoeff binomial_coeff;
        binomial_coeff.construct(p);
        for (int32 ii = 0; ii <= p; ii++)
        {
          m_table[ii] = binomial_coeff.get_val();
          binomial_coeff.slide_over(0);
        }
      }

      DRAY_EXEC const int32 & operator[](int32 i) const { return m_table[i]; }
    };

    /** eval_d_edge(ShapeHex, General) */
    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d_edge( ShapeHex,
                                             const OrderPolicy<General> order_p,
                                             const int32 eid,
                                             const ReadDofPtr<Vec<Float, ncomp>> &C,
                                             const Vec<Float, 1> &rc,
                                             Vec<Vec<Float, ncomp>, 1> &out_deriv )
    {
      const int32 p = eattr::get_order(order_p);
      const HexEdgeWalker<General> hew(order_p, eid);
      const Float &u = rc[0];
      const Float ubar = 1.0 - u;

      if (p == 0)
      {
        out_deriv[0] = 0;
        return C[0];
      }

      BinomialCoeffTable B(p-1);

      Float upow = 1.0;
      Vec<Float, ncomp> val_u_L, val_u_R;
      val_u_L = 0;
      val_u_R = 0;

      int32 i = 0;
      Vec<Float, ncomp> Ci = C[ hew.m_edge_base + 0 * hew.m_edge_stride ];
      for (i = 1; i <= p; ++i)
      {
        val_u_L = val_u_L * ubar + Ci * (B[i-1] * upow);
        Ci = C[ hew.m_edge_base + i * hew.m_edge_stride ];
        val_u_R = val_u_R * ubar + Ci * (B[i-1] * upow);
        upow *= u;
      }

      out_deriv[0] = (val_u_R - val_u_L) * p;
      return (val_u_R * u) + (val_u_L * ubar);
    }


    /** eval_d_face(ShapeHex, Linear) */
    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d_face( ShapeHex,
                                             const OrderPolicy<Linear> order_p,
                                             const int32 fid,
                                             const ReadDofPtr<Vec<Float, ncomp>> &C,
                                             const Vec<Float, 2> &rc,
                                             Vec<Vec<Float, ncomp>, 2> &out_deriv )
    {
      constexpr int32 p = eattr::get_order(order_p.as_cxp());
      const HexFaceWalker<Linear> hfw(order_p.as_cxp(), fid);

      const Float &u = rc[0],  _u = 1.0-u;
      const Float &v = rc[1],  _v = 1.0-v;

      const Vec<Float, ncomp> C00 = C[ hfw.m_f_base + 0 * hfw.m_f_stride1 + 0 * hfw.m_f_stride0 ];
      const Vec<Float, ncomp> C01 = C[ hfw.m_f_base + 0 * hfw.m_f_stride1 + 1 * hfw.m_f_stride0 ];
      const Vec<Float, ncomp> C10 = C[ hfw.m_f_base + 1 * hfw.m_f_stride1 + 0 * hfw.m_f_stride0 ];
      const Vec<Float, ncomp> C11 = C[ hfw.m_f_base + 1 * hfw.m_f_stride1 + 1 * hfw.m_f_stride0 ];

      out_deriv[0] = ((C11-C10)*p)*v + ((C01-C00)*p)*_v;
      out_deriv[1] = ((C11-C01)*p)*u + ((C10-C00)*p)*_u;
      return (C11*u + C10*_u)*v + (C01*u + C00*_u)*_v;
    }

    /** eval_d_face(ShapeHex, Quadratic) */
    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d_face( ShapeHex,
                                             const OrderPolicy<Quadratic> order_p,
                                             const int32 fid,
                                             const ReadDofPtr<Vec<Float, ncomp>> &C,
                                             const Vec<Float, 2> &rc,
                                             Vec<Vec<Float, ncomp>, 2> &out_deriv )
    {
      constexpr int32 p = eattr::get_order(order_p.as_cxp());
      const HexFaceWalker<Quadratic> hfw(order_p.as_cxp(), fid);

      const Float &u = rc[0],  _u = 1.0-u;
      const Float &v = rc[1],  _v = 1.0-v;

      Vec<Float, ncomp> C00 = C[ hfw.m_f_base + 0 * hfw.m_f_stride1 + 0 * hfw.m_f_stride0 ];
      Vec<Float, ncomp> C01 = C[ hfw.m_f_base + 0 * hfw.m_f_stride1 + 1 * hfw.m_f_stride0 ];
      Vec<Float, ncomp> C02 = C[ hfw.m_f_base + 0 * hfw.m_f_stride1 + 2 * hfw.m_f_stride0 ];
      C00 = C01*u + C00*_u;  //DeCasteljau
      C01 = C02*u + C01*_u;  //DeCasteljau

      Vec<Float, ncomp> C10 = C[ hfw.m_f_base + 1 * hfw.m_f_stride1 + 0 * hfw.m_f_stride0 ];
      Vec<Float, ncomp> C11 = C[ hfw.m_f_base + 1 * hfw.m_f_stride1 + 1 * hfw.m_f_stride0 ];
      Vec<Float, ncomp> C12 = C[ hfw.m_f_base + 1 * hfw.m_f_stride1 + 2 * hfw.m_f_stride0 ];
      C10 = C11*u + C10*_u;  //DeCasteljau
      C11 = C12*u + C11*_u;  //DeCasteljau

      Vec<Float, ncomp> C20 = C[ hfw.m_f_base + 2 * hfw.m_f_stride1 + 0 * hfw.m_f_stride0 ];
      Vec<Float, ncomp> C21 = C[ hfw.m_f_base + 2 * hfw.m_f_stride1 + 1 * hfw.m_f_stride0 ];
      Vec<Float, ncomp> C22 = C[ hfw.m_f_base + 2 * hfw.m_f_stride1 + 2 * hfw.m_f_stride0 ];
      C20 = C21*u + C20*_u;  //DeCasteljau
      C21 = C22*u + C21*_u;  //DeCasteljau

      C00 = C10*v + C00*_v;  //DeCasteljau
      C10 = C20*v + C10*_v;  //DeCasteljau

      C01 = C11*v + C01*_v;  //DeCasteljau
      C11 = C21*v + C11*_v;  //DeCasteljau

      out_deriv[0] = ((C11-C10)*p)*v + ((C01-C00)*p)*_v;
      out_deriv[1] = ((C11-C01)*p)*u + ((C10-C00)*p)*_u;
      return (C11*u + C10*_u)*v + (C01*u + C00*_u)*_v;
    }

    namespace detail
    {
      DRAY_EXEC Float shape(const BinomialCoeffTable &B,
                            int32 p,
                            int32 i,
                            const Float &u,
                            const Float &_u)
      {
        return B[i] * ipow_w(_u, p-i) * ipow_w(u, i);
      }

      DRAY_EXEC Float shape(const BinomialCoeffTable &B,
                            int32 p,
                            int32 i,
                            const Float &u,
                            const Float &_u,
                            Float &dshape)
      {
        dshape = B[i] * ( (i==p ? 0 : -(p-i) * ipow_w(_u, p-i-1) * ipow_w(u, i))
                        + (i==0 ? 0 :      i * ipow_w(_u, p-i)   * ipow_w(u, i-1)) );

        return B[i] * ipow_w(_u, p-i) * ipow_w(u, i);
      }
    }

    /** eval_d_face(ShapeHex, General) */
    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d_face( ShapeHex,
                                             const OrderPolicy<General> order_p,
                                             const int32 fid,
                                             const ReadDofPtr<Vec<Float, ncomp>> &C,
                                             const Vec<Float, 2> &rc,
                                             Vec<Vec<Float, ncomp>, 2> &out_deriv )
    {
      const int32 p = eattr::get_order(order_p);
      const HexFaceWalker<General> hfw(order_p, fid);
      const ReadDofPtr<Vec<Float, ncomp>> CF = C + hfw.m_f_base;
      const int32 &stride0 = hfw.m_f_stride0;
      const int32 &stride1 = hfw.m_f_stride1;

      if (p == 0)
      {
        out_deriv[0] = 0;
        return C[0];
      }

      BinomialCoeffTable B(p);

      const Float &u = rc[0],  _u = 1.0-u;
      const Float &v = rc[1],  _v = 1.0-v;

      Vec<Float, ncomp> result;
      result = 0;
      out_deriv = 0;

      // Power rule version -- simple.
      //   (Horner's rule version is just too hard to read.)
      for (int32 j = 0; j <= p; ++j)
      {
        Float shape_j, dshape_j;
        shape_j = detail::shape(B, p, j, v, _v, dshape_j);

        for (int32 i = 0; i <= p; ++i)
        {
          Float shape_i, dshape_i;
          shape_i = detail::shape(B, p, i, u, _u, dshape_i);

          const Vec<Float, ncomp> C_ij = CF[ j * stride1 + i * stride0 ];

          result += C_ij * shape_i * shape_j;
          out_deriv[0] += C_ij * dshape_i * shape_j;
          out_deriv[1] += C_ij * shape_i * dshape_j;
        }
      }

      return result;
    }






    /** eval_d() */
    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d( ShapeTri,
                                        OrderPolicy<Linear>,
                                        const ReadDofPtr<Vec<Float, ncomp>> &C,
                                        const Vec<Float, 2> &rc,
                                        Vec<Vec<Float, ncomp>, 2> &out_deriv )
    {
      // C[2]
      // C[0] C[1]
      const Float &u = rc[0], &v = rc[1];
      const Float t = 1.0f - u - v;

      Float sd[3];
      sd[0] = -1.0f;

      // du
      sd[1] = 1.0f;   sd[2] = 0.0f;
      out_deriv[0] = C[0] * sd[0] + C[1] * sd[1] + C[2] * sd[2];

      // du
      sd[1] = 0.0f;   sd[2] = 1.0f;
      out_deriv[0] = C[0] * sd[0] + C[1] * sd[1] + C[2] * sd[2];

      const Float s[3] = { t, u, v };
      return C[0] * s[0] + C[1] * s[1] + C[2] * s[2];
    }

    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d( ShapeTri,
                                        OrderPolicy<Quadratic>,
                                        const ReadDofPtr<Vec<Float, ncomp>> &C,
                                        const Vec<Float, 2> &rc,
                                        Vec<Vec<Float, ncomp>, 2> &out_deriv )
    {
      // C[6]
      //
      // C[3] C[4]
      //
      // C[0] C[1] C[2]
      const Float &u = rc[0], &v = rc[1];
      const Float t = 1.0f - u - v;

      Float sd[6];
      sd[0] = 2*(-t);

      // -------------------------------
      // du
                        sd[1] = 2*(t-u);    sd[2] = 2*u;
      sd[3] = 2*(-v);   sd[4] = 2*(v);
      sd[5] = 0.0f;
      //
      out_deriv[0] = C[0] * sd[0] + C[1] * sd[1] + C[2] * sd[2]
                   + C[3] * sd[3] + C[4] * sd[4]
                   + C[5] * sd[5];
      // -------------------------------

      // -------------------------------
      // dv
                         sd[1] = 2*(-u);   sd[2] = 0.0f;
      sd[3] = 2*(t-v);   sd[4] = 2*(u);
      sd[5] = 2*v;
      //
      out_deriv[1] = C[0] * sd[0] + C[1] * sd[1] + C[2] * sd[2]
                   + C[3] * sd[3] + C[4] * sd[4]
                   + C[5] * sd[5];
      // -------------------------------


      const Float s[6] = { t*t,     2*t*u,   u*u,
                           2*t*v,   2*v*u,
                           v*v };
      return C[0] * s[0] + C[1] * s[1] + C[2] * s[2] +
             C[3] * s[3] + C[4] * s[4] +
             C[5] * s[5];
    }



    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d( ShapeTet,
                                        OrderPolicy<Linear>,
                                        const ReadDofPtr<Vec<Float, ncomp>> &C,
                                        const Vec<Float, 3> &rc,
                                        Vec<Vec<Float, ncomp>, 3> &out_deriv )
    {
      //  C[2]
      //
      //  C[0]  C[1]
      // C[3]
      const Float &u = rc[0], &v = rc[1], &w = rc[2];
      const Float t = 1.0f - u - v - w;

      Float sd[4];
      sd[0] = -1.0f;

      // du
      sd[1] = 1.0f;   sd[2] = 0.0f;   sd[3] = 0.0f;
      out_deriv[0] = C[0] * sd[0] + C[1] * sd[1] + C[2] * sd[2];

      // du
      sd[1] = 0.0f;   sd[2] = 1.0f;   sd[3] = 0.0f;
      out_deriv[0] = C[0] * sd[0] + C[1] * sd[1] + C[2] * sd[2];

      // dw
      sd[1] = 0.0f;   sd[2] = 0.0f;   sd[3] = 1.0f;
      out_deriv[0] = C[0] * sd[0] + C[1] * sd[1] + C[2] * sd[2];

      const Float s[4] = { t, u, v, w };
      return C[0] * s[0] + C[1] * s[1] + C[2] * s[2] + C[3] * s[3];
    }



    template <int32 ncomp>
    DRAY_EXEC Vec<Float, ncomp> eval_d( ShapeTet,
                                        OrderPolicy<Quadratic>,
                                        const ReadDofPtr<Vec<Float, ncomp>> &C,
                                        const Vec<Float, 3> &rc,
                                        Vec<Vec<Float, ncomp>, 3> &out_deriv )
    {
      //  Behold, the P2 tetrahedron
      //
      //              v=1
      //
      //              C[5]
      //             /|   `
      //            / C[3]  C[4]
      //           C[8]        `
      //          /|  C[0]--C[1]--C[2]   u=1
      //         / C[6]   C[7]  '
      //        C[9]   '
      //    w=1
      //
      const Float &u = rc[0], &v = rc[1], &w = rc[2];
      const Float t = 1.0f - u - v - w;

      Float sd[10];
      sd[0] = 2*(-t);

      // -------------------------------
      // du
                        sd[1] = 2*(t-u);   sd[2] = 2*u;
      sd[3] = 2*(-v);   sd[4] = 2*(v);
      sd[5] = 0.0f;
                          sd[6] = 2*(-w);    sd[7] = 2*(w);
                          sd[8] = 0.0f;
                                                sd[9] = 0.0f;
      out_deriv[0] = 0;
      for (int32  i = 0; i < 10; ++i)
        out_deriv[0] += C[i] * sd[i];
      // -------------------------------

      // -------------------------------
      // dv
                         sd[1] = 2*(-u);  sd[2] = 0.0f;
      sd[3] = 2*(t-v);   sd[4] = 2*(u);
      sd[5] = 2*v;
                          sd[6] = 2*(-w);    sd[7] = 0.0f;
                          sd[8] = 2*(w);
                                                sd[9] = 0.0f;
      out_deriv[1] = 0;
      for (int32  i = 0; i < 10; ++i)
        out_deriv[1] += C[i] * sd[i];
      // -------------------------------

      // -------------------------------
      // dv

                       sd[1] = 2*(-u);   sd[2] = 0.0f;
      sd[3] = 2*(-v);  sd[4] = 0.0f;
      sd[5] = 0.0f;
                         sd[6] = 2*(t-w);   sd[7] = 2*(u);
                         sd[8] = 2*(v);
                                              sd[9] = 2*w;
      out_deriv[2] = 0;
      for (int32  i = 0; i < 10; ++i)
        out_deriv[2] += C[i] * sd[i];
      // -------------------------------


      const Float s[10] = { t*t,     2*t*u,   u*u,
                            2*t*v,   2*v*u,
                            v*v,
                                       2*t*w,  2*u*w,
                                       2*v*w,
                                                  w*w };

      Vec<Float, ncomp> ret;  ret = 0;
      for (int32 i = 0; i < 10; ++i)
        ret += C[i] * s[i];
      return ret;
    }


  }//eops



  // --------------------------------------------------------------------------
  // split_inplace()
  // --------------------------------------------------------------------------

  namespace detail
  {
    constexpr int32 cartesian_to_tri_idx(int32 i, int32 j, int32 edge);
    constexpr int32 cartesian_to_tet_idx(int32 i, int32 j, int32 k, int32 e);
  }

  // The Split<> object describes a binary split of the simplex at some point
  // (given by 'factor') along an edge (given by vtx_displaced and vtx_tradeoff).
  // Each row of coefficients parallel to the specified edge undergoes
  // 1D-DeCasteljau subdivision. The side closest to the 'tradeoff' vertex is
  // treated as the fixed, exterior node. The side closest to the 'displaced'
  // vertex is treated as the parameter-dependent, interior node.
  //
  //              .                 .           .           .
  //             .-*               . .         . .         . .
  //            .-*-*             . .-*       . . .       . . .
  //           .-*-*-*           . .-*-*     . . .-*     . . . .
  //       (v1=p)    (v0=p)
  //     tradeoff    displaced
  //
  // Subject to axis permutations, the splitting can be carried out as below:
  //
  //   // Triangle
  //   for (0 <= h < p)
  //     for (p-h >= num_updates >= 1)
  //       for (p-h >= v0 > p-h - num_updates, v0+v1 = p-h)
  //         C[v0,v1;h] := f*C[v0,v1;h] + (1-f)*C[v0-1,v1+1;h];
  //
  //   // Tetrahedron
  //   for (0 <= g+h < p)
  //     for (p-(g+h) >= num_updates >= 1)
  //       for (p-(g+h) >= v0 > p-(g+h) - num_updates, v0+v1 = p-(g+h))
  //         C[v0,v1;g,h] := f*C[v0,v1;g,h] + (1-f)*C[v0-1,v1+1;g,h];
  //

  //
  // split_inplace<2, Simplex>        (Triangle)
  //
  template <int32 ncomp, int32 P>
  DRAY_EXEC void split_inplace(
      ShapeTri, OrderPolicy<P> order_p,
      WriteDofPtr<Vec<Float, ncomp>> dof_ptr,
      const Split<ElemType::Simplex> &split)
  {
    const uint8 p = (uint8) eattr::get_order(order_p);
    const uint8 v0 = (uint8) split.vtx_displaced;
    const uint8 v1 = (uint8) split.vtx_tradeoff;
    const uint8 v2 = 0+1+2 - v0 - v1;

    uint8 b[3];  // barycentric indexing

    // I think this way of expressing the permuation is most readable.
    // On the other hand, potential for loop unrolling might be easier to
    // detect if the permutation was expressed using the inverse.

    for (b[v2] = 0; b[v2] < p; ++b[v2])
      for (uint8 num_updates = p-b[v2]; num_updates >= 1; --num_updates)
        for (b[v0] = p-b[v2]; b[v0] > p-b[v2] - num_updates; --b[v0])
        {
          b[v1] = p-b[v2]-b[v0];

          uint8 b_left[3];
          b_left[v0] = b[v0] - 1;
          b_left[v1] = b[v1] + 1;
          b_left[v2] = b[v2];

          const uint32 right = detail::cartesian_to_tri_idx(b[0], b[1], p+1);
          const uint32 left = detail::cartesian_to_tri_idx(b_left[0], b_left[1], p+1);

          dof_ptr[right] =
              dof_ptr[right] * split.factor + dof_ptr[left] * (1-split.factor);
        }
  }

  /** @deprecated */
  template <int32 ncomp, int32 P>
  DRAY_EXEC void split_inplace(
      const Element<2, ncomp, ElemType::Simplex, P> &elem_info,  // tag for template + order
      WriteDofPtr<Vec<Float, ncomp>> dof_ptr,
      const Split<ElemType::Simplex> &split)
  {
    split_inplace(ShapeTri{},
                  eattr::adapt_create_order_policy(OrderPolicy<P>{}, elem_info.get_order()),
                  dof_ptr,
                  split);
  }



  //
  // split_inplace<3, Simplex>          (Tetrahedron)
  //
  template <int32 ncomp, int32 P>
  DRAY_EXEC void split_inplace(
      ShapeTet, OrderPolicy<P> order_p,
      WriteDofPtr<Vec<Float, ncomp>> dof_ptr,
      const Split<ElemType::Simplex> &split)
  {
    const uint8 p = (uint8) eattr::get_order(order_p);
    const uint8 v0 = (uint8) split.vtx_displaced;
    const uint8 v1 = (uint8) split.vtx_tradeoff;

    const uint8 avail = -1u & ~(1u << v0) & ~(1u << v1);
    const uint8 v2 = (avail & 1u) ? 0 : (avail & 2u) ? 1 : (avail & 4u) ? 2 : 3;
    const uint8 v3 = 0+1+2+3 - v0 - v1 - v2;

    uint8 b[4];  // barycentric indexing

    for (b[v3] = 0; b[v3] < p; ++b[v3])
      for (b[v2] = 0; b[v2] < p-b[v3]; ++b[v2])
      {
        const uint8 gph = b[v2] + b[v3];

        for (uint8 num_updates = p-gph; num_updates >= 1; --num_updates)
          for (b[v0] = p-gph; b[v0] > p-gph - num_updates; --b[v0])
          {
            b[v1] = p-gph-b[v0];

            int8 b_left[4];
            b_left[v0] = b[v0] - 1;
            b_left[v1] = b[v1] + 1;
            b_left[v2] = b[v2];
            b_left[v3] = b[v3];

            const uint32 right = detail::cartesian_to_tet_idx(
                b[0], b[1], b[2], p+1);
            const uint32 left = detail::cartesian_to_tet_idx(
                b_left[0], b_left[1], b_left[2], p+1);

            dof_ptr[right] =
                dof_ptr[right] * split.factor + dof_ptr[left] * (1-split.factor);
          }
      }
  }

  /** @deprecated */
  template <int32 ncomp, int32 P>
  DRAY_EXEC void split_inplace(
      const Element<3, ncomp, ElemType::Simplex, P> &elem_info,  // tag for template + order
      WriteDofPtr<Vec<Float, ncomp>> dof_ptr,
      const Split<ElemType::Simplex> &split)
  {
    split_inplace(ShapeTet{},
                  eattr::adapt_create_order_policy(OrderPolicy<P>{}, elem_info.get_order()),
                  dof_ptr,
                  split);
  }

  // Binary split on quad:
  //
  //  left:
  //     .-*-*-*    . .-*-*    . . .-*    . . . .
  //     .-*-*-*    . .-*-*    . . .-*    . . . .
  //     .-*-*-*    . .-*-*    . . .-*    . . . .
  //     .-*-*-*    . .-*-*    . . .-*    . . . .
  //
  //  right:
  //     *-*-*-.    *-*-. .    *-. . .    . . . .
  //     *-*-*-.    *-*-. .    *-. . .    . . . .
  //     *-*-*-.    *-*-. .    *-. . .    . . . .
  //     *-*-*-.    *-*-. .    *-. . .    . . . .
  //

  //
  // split_inplace<Tensor>
  //
  template <int32 ncomp, int32 P>
  DRAY_EXEC void split_inplace(
      ShapeHex, OrderPolicy<P> order_p,
      WriteDofPtr<Vec<Float, ncomp>> dof_ptr,
      const Split<ElemType::Tensor> &split)
  {
    constexpr int32 dim = eattr::get_dim(ShapeHex{});
    const uint32 p = eattr::get_order(order_p);

    uint32 p_order_pow[4];
    p_order_pow[0] = 1;
    p_order_pow[1] = p_order_pow[0] * (p + 1);
    p_order_pow[2] = p_order_pow[1] * (p + 1);
    p_order_pow[3] = p_order_pow[2] * (p + 1);

    const int32 &axis = split.axis;
    static_assert((1 <= dim && dim <= 3), "split_inplace() only supports 1D, 2D, or 3D.");
    assert((0 <= axis && axis < dim));
    const uint32 stride = p_order_pow[axis];
    const uint32 chunk_sz = p_order_pow[axis+1];
    const uint32 num_chunks = p_order_pow[dim - (axis+1)];

    if (!split.f_lower_t_upper)
    {
      // Left
      for (int32 chunk = 0; chunk < num_chunks; ++chunk, dof_ptr += chunk_sz)
      {
        // Split the chunk along axis.
        // If there are axes below axis, treat them as a vector of dofs.

        // In DeCasteljau left split, we repeatedly overwrite the right side.
        for (int32 from_front = 1; from_front <= p; ++from_front)
          for (int32 ii = p; ii >= 0+from_front; --ii)
            for (int32 e = 0; e < stride; ++e)
            {
              dof_ptr[ii*stride + e] =
                  dof_ptr[(ii-1)*stride + e] * (1 - split.factor)
                  + dof_ptr[ii*stride + e] * (split.factor);
            }
      }
    }
    else
    {
      // Right
      for (int32 chunk = 0; chunk < num_chunks; ++chunk, dof_ptr += chunk_sz)
      {
        // Split the chunk along axis.
        // If there are axes below axis, treat them as a vector of dofs.

        // In DeCasteljau right split, we repeatedly overwrite the left side.
        for (int32 from_back = 1; from_back <= p; ++from_back)
          for (int32 ii = 0; ii <= p-from_back; ++ii)
            for (int32 e = 0; e < stride; ++e)
            {
              dof_ptr[ii*stride + e] =
                  dof_ptr[ii*stride + e] * (1 - split.factor)
                  + dof_ptr[(ii+1)*stride + e] * (split.factor);
            }
      }
    }
  }

  //
  // split_inplace<Tensor>
  //
  template <int32 ncomp, int32 P>
  DRAY_EXEC void split_inplace(
      ShapeQuad, OrderPolicy<P> order_p,
      WriteDofPtr<Vec<Float, ncomp>> dof_ptr,
      const Split<ElemType::Tensor> &split)
  {
    constexpr int32 dim = eattr::get_dim(ShapeQuad{});
    const uint32 p = eattr::get_order(order_p);

    uint32 p_order_pow[4];
    p_order_pow[0] = 1;
    p_order_pow[1] = p_order_pow[0] * (p + 1);
    p_order_pow[2] = p_order_pow[1] * (p + 1);
    p_order_pow[3] = p_order_pow[2] * (p + 1);

    const int32 &axis = split.axis;
    static_assert((1 <= dim && dim <= 3), "split_inplace() only supports 1D, 2D, or 3D.");
    assert((0 <= axis && axis < dim));
    const uint32 stride = p_order_pow[axis];
    const uint32 chunk_sz = p_order_pow[axis+1];
    const uint32 num_chunks = p_order_pow[dim - (axis+1)];

    if (!split.f_lower_t_upper)
    {
      // Left
      for (int32 chunk = 0; chunk < num_chunks; ++chunk, dof_ptr += chunk_sz)
      {
        // Split the chunk along axis.
        // If there are axes below axis, treat them as a vector of dofs.

        // In DeCasteljau left split, we repeatedly overwrite the right side.
        for (int32 from_front = 1; from_front <= p; ++from_front)
          for (int32 ii = p; ii >= 0+from_front; --ii)
            for (int32 e = 0; e < stride; ++e)
            {
              dof_ptr[ii*stride + e] =
                  dof_ptr[(ii-1)*stride + e] * (1 - split.factor)
                  + dof_ptr[ii*stride + e] * (split.factor);
            }
      }
    }
    else
    {
      // Right
      for (int32 chunk = 0; chunk < num_chunks; ++chunk, dof_ptr += chunk_sz)
      {
        // Split the chunk along axis.
        // If there are axes below axis, treat them as a vector of dofs.

        // In DeCasteljau right split, we repeatedly overwrite the left side.
        for (int32 from_back = 1; from_back <= p; ++from_back)
          for (int32 ii = 0; ii <= p-from_back; ++ii)
            for (int32 e = 0; e < stride; ++e)
            {
              dof_ptr[ii*stride + e] =
                  dof_ptr[ii*stride + e] * (1 - split.factor)
                  + dof_ptr[(ii+1)*stride + e] * (split.factor);
            }
      }
    }
  }

  /** @deprecated */
  template <int32 dim, int32 ncomp, int32 P>
  DRAY_EXEC void split_inplace(
      const Element<dim, ncomp, ElemType::Tensor, P> &elem_info,  // tag for template + order
      WriteDofPtr<Vec<Float, ncomp>> dof_ptr,
      const Split<ElemType::Tensor> &split)
  {
    split_inplace(Shape<dim, Tensor>{},
                  eattr::adapt_create_order_policy(OrderPolicy<P>{}, elem_info.get_order()),
                  dof_ptr,
                  split);

  }

  // --------------------------------------------------------------------------



  // Imported from isosurface_meshing

/// namespace DeCasteljauSplitting
/// {
/// // sub_element()
/// //
/// // Replaces sub_element_fixed_order()
/// // This version operates in-place and does not assume fixed order.
/// template <uint32 dim, uint32 ncomp>
/// DRAY_EXEC void sub_element(uint32 p_order,
///                            const Range *ref_boxs,
///                            WriteDofPtr<Vec<Float, ncomp>> &wptr)
/// {
///   // Split along each axis sequentially. It is a tensor element.
///   for (int32 d = 0; d < dim; ++d)
///   {
///     const Float t1 = ref_boxs[d].max();
///     Float t0 = ref_boxs[d].min();
///
///     // Split left, using right endpoint.
///     if (t1 < 1.0)
///       split_inplace_left(wptr, t1, dim, d, p_order);
///
///     // Left endpoint relative to the new subinterval.
///     if (t1 > 0.0) t0 /= t1;
///
///     // Split right, using left endpoint.
///     if (t0 > 0.0)
///       split_inplace_right(wptr, t0, dim, d, p_order);
///   }
/// }





}//namespace dray

#endif//DRAY_ELEM_OPS_HPP
