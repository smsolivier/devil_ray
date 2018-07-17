#include "gtest/gtest.h"
#include <dray/high_order_shape.hpp>
#include <dray/array.hpp>
#include <dray/vec.hpp>
#include <dray/matrix.hpp>
#include <dray/math.hpp>
#include <dray/binomial.hpp>
#include <dray/utils/png_encoder.hpp>

#include <iostream>
#include <string>
#include <string.h>


TEST(dray_test, dray_high_order_shape)
{

  // -- Test SimpleTensor -- //
  {
    constexpr int t_order = 3;
    constexpr int vec_size = 3;
    dray::SimpleTensor<float> st;
    st.s = vec_size;
    const int size_tensor = st.get_size_tensor(t_order);
    float *some_mem = new float[size_tensor];
    float * ptrs3[t_order];
    
    // Store vector data.
    st.get_vec_init_ptrs(t_order, some_mem, ptrs3);
    ptrs3[0][0] = 1;
    ptrs3[1][0] = 5;
    ptrs3[2][0] = 8;
    for (int ii = 1; ii < vec_size; ii++)
    {
      ptrs3[0][ii] = ptrs3[0][ii-1] + 1;
      ptrs3[1][ii] = ptrs3[1][ii-1] + 1;
      ptrs3[2][ii] = ptrs3[2][ii-1] + 1;
    }

    // Construct.
    st.construct_in_place(t_order, some_mem);

    // Output results.
    std::cout << "The tensor is..." << std::endl;
    for (int ii = 0; ii < size_tensor; ii++)
    {
      printf("%.f  ", some_mem[ii]);
    }
    printf("\n");

    delete [] some_mem;
  }


  // -- Test PowerBasis and BernsteinBasis -- //
  // 1D cubic.
  {
    std::cout << "Cubic Univariate." << std::endl;

    constexpr int power = 3;
    float control_vals[power+1] = {22, 10, 5, 2};
    
    constexpr int num_samples = 4;
    ////float sample_vals[num_samples] = {0, 1, 2, 3};
    float sample_vals[num_samples] = {0, .3, .6, .9};

    // Print symbolic polynomial.
    printf("P(x) = (%f)", control_vals[0]);
    for(int k = 1; k <= power; k++)
      printf(" + (%f)x^%d", control_vals[k], k);
    printf("\n");

    dray::Vec<float,1> x;
    typedef dray::Vec<float,1> CtrlPtType;
    const CtrlPtType *coeff = (const CtrlPtType *) control_vals;
    dray::Vec<float,1> val;
    dray::Vec<dray::Vec<float,1>,1> deriv;
    for (int ii = 0; ii < num_samples; ii++)
    {
      x = sample_vals[ii];
      dray::PowerBasis<float,1> pb;
      pb.init_shape(power);
      pb.linear_combo<const CtrlPtType*, 1>(x, coeff, val, deriv);
      printf("P(%f)  = %f\n", x[0], val[0]);
      printf("P'(%f) = %f\n", x[0], deriv[0][0]);
    }

    for (int ii = 0; ii < num_samples; ii++)
    {
      x = sample_vals[ii];
      float aux_mem[8];
      dray::BernsteinBasis<float,1> bb(power, aux_mem);
      printf("bb.get_aux_req() == %d\n", bb.get_aux_req());
      bb.linear_combo<const CtrlPtType*, 1>(x, coeff, val, deriv);
      std::cout << "Bernstein basis:    " << val << deriv << std::endl;
    }
  }

  // 3D linear.
  {
    std::cout << "Linear Trivariate." << std::endl;

    constexpr int power = 1;

    float control_vals[8] = {1,2,3,4,5,6,7,8};

    //dray::Vec<float,3> xyz = {1,2,3};
    dray::Vec<float,3> xyz = {.5,.5,.5};
    typedef dray::Vec<float,1> CtrlPtType;
    const CtrlPtType *coeff = (const CtrlPtType *) control_vals;
    dray::Vec<float,1> val;
    //dray::Matrix<float,1,3> deriv;
    dray::Vec<dray::Vec<float,1>,3> deriv;

    dray::PowerBasis<float,3> pb;
    pb.init_shape(power);
    pb.linear_combo<const CtrlPtType*, 1>(xyz, coeff, val, deriv);
    std::cout << "P(" << xyz << ") = " << val[0] << std::endl;
    std::cout << "P'(" << xyz << ") = " << deriv << std::endl;

    float aux_mem[12];
    dray::BernsteinBasis<float,3> bb(power, aux_mem);
    printf("bb.get_aux_req() == %d\n", bb.get_aux_req());
    bb.linear_combo<const CtrlPtType*, 1>(xyz, coeff, val, deriv);
    std::cout << "Bernstein basis:    " << val << deriv << std::endl;
  }

  // 3D quadratic.
  {
    std::cout << "Quadratic Trivariate." << std::endl;

    constexpr int power = 2;

    float control_vals[27] = {1, 2, 3, 4, 5, 6, 7, 8, 9,
                              10,11,12,13,14,15,16,17,18,
                              19,20,21,22,23,24,25,26,27};

    /////dray::Vec<float,3> xyz = {1, 2, 3};
    dray::Vec<float,3> xyz = {.33, .5, .66};
    ////dray::Vec<float,3> xyz = {.5,.5,.5};
    typedef dray::Vec<float,1> CtrlPtType;
    const CtrlPtType *coeff = (const CtrlPtType *) control_vals;
    dray::Vec<float,1> val;
    dray::Vec<dray::Vec<float,1>,3> deriv;

    dray::PowerBasis<float,3> pb;
    pb.init_shape(power);
    pb.linear_combo<const CtrlPtType*, 1>(xyz, coeff, val, deriv);
    std::cout << "P(" << xyz << ") = " << val[0] << std::endl;
    std::cout << "P'(" << xyz << ") = " << deriv << std::endl;

    // Check
    {
      float x = xyz[0], y = xyz[1], z = xyz[2];
      int coeff_idx = 0;
      float xa = 0;
      for (int xi = power; xi >= 0; xi--)
      {
        float ya = 0;
        int cx = coeff_idx + 9*xi;
        for (int yi = power; yi >= 0; yi--)
        {
          float za = 0;
          int cy = cx + 3*yi;
          for (int zi = power; zi >= 0; zi--)
          {
            int cz = cy + zi;
            za = za * z + control_vals[cz];
          }  // zi
          ya = ya * y + za;
        }  // yi
        xa = xa * x + ya;
      }  // xi
      std::cout << "P(x,y,z) should be    " << xa << std::endl;
    }

    float aux_mem[18];
    dray::BernsteinBasis<float,3> bb(power, aux_mem);
    printf("bb.get_aux_req() == %d\n", bb.get_aux_req());
    bb.linear_combo<const CtrlPtType*, 1>(xyz, coeff, val, deriv);
    std::cout << "Bernstein basis:    " << val << deriv << std::endl;

    // Check by direct computation (recomputing powers many times, fixed at 3 dimensions).
    {
      float binom[power+1];
      dray::BinomRow<float>::fill_single_row(power, binom);
      float x = xyz[0], y = xyz[1], z = xyz[2];
      float val = 0;
      float deriv[3] = {0,0,0};
      int dof_idx = 0;
      for (int ii = 0; ii <= power; ii++)
      {
        float x_shape = binom[ii] * pow(x,ii) * pow(1.0-x, power-ii);
        float dx_shape = binom[ii] * ( ii==0 ? -power * pow(1.0-x,power-1) :
                                       ii==power ? power * pow(x,power-1) :
                                           (ii-power*x)*pow(x,ii-1)*pow(1.-x,power-ii-1) );
        for (int jj = 0; jj <= power; jj++)
        {
          float y_shape = binom[jj] * pow(y,jj) * pow(1.0-y, power-jj);
          float dy_shape = binom[jj] * ( jj==0 ? -power * pow(1.0-y,power-1) :
                                         jj==power ? power * pow(y,power-1) :
                                             (jj-power*y)*pow(y,jj-1)*pow(1.-y,power-jj-1) );
          for (int kk = 0; kk <= power; kk++)
          {
            float z_shape = binom[kk] * pow(z,kk) * pow(1.0-z, power-kk);
            float dz_shape = binom[kk] * ( kk==0 ? -power * pow(1.0-z,power-1) :
                                           kk==power ? power * pow(z,power-1) :
                                               (kk-power*z)*pow(z,kk-1)*pow(1.-z,power-kk-1) );
            
            float ctrl_val = control_vals[dof_idx++];
            val += ctrl_val * x_shape * y_shape * z_shape;
            deriv[0] += ctrl_val * dx_shape * y_shape * z_shape;
            deriv[1] += ctrl_val * x_shape * dy_shape * z_shape;
            deriv[2] += ctrl_val * x_shape * y_shape * dz_shape;
          }
        }
      }
      printf("Should be val == %f    deriv == (%f  %f  %f)\n", val, deriv[0],deriv[1],deriv[2]);
    }
  }

  //--- Test QuerySum and QueryCat---//

  //--- Test linear shape evaluator---//    Linear Shape works.
  //--- Test Bernstein shape evaluator---//
  //--- Test ElTrans ---//
   
    // There are two quadratic unit-cubes, adjacent along X, sharing a face in the YZ plane.
    // There are 45 total control points: 2 vol mids, 11 face mids, 20 edge mids, and 12 vertices.
    float grid_vals[45] = 
        { 10, -10,                           // 0..1 vol mids A and B
          15,7,7,7,7,  0, -15,-7,-7,-7,-7,      // 2..12 face mids A(+X,+Y,+Z,-Y,-Z) AB B(-X,+Y,+Z,-Y,-Z)
          12,12,12,12,  -12,-12,-12,-12,     // 13..20 edge mids on ends +X/-X A(+Y,+Z,-Y,-Z) B(+Y,+Z,-Y,-Z)
          5,5,5,5,  -5,-5,-5,-5,             // 21..28 edge mids YZ corners  A(++,-+,--,+-) B(++,-+,--,+-)
          0,0,0,0,                           // 29..32 edge mids on shared face AB(+Y,+Z,-Y,-Z)
          20,20,20,20,  -20,-20,-20,-20,     // 33..40 vertices on ends +X/-X, YZ corners A(++,-+,--,+-) B(++,-+,--,+-)
          0,0,0,0 };                         // 41..44 vertices on shared face, YZ corners AB(++,-+,--,+-)

    // Map the per-element degrees of freedom into the total set of control points.
    int ctrl_idx[54];
    int * const ax = ctrl_idx, * const bx = ctrl_idx + 27;

    // Nonshared nodes.
    ax[13] = 0;  bx[13] = 1;

    ax[22] = 2;  bx[4] = 8;
    ax[16] = 3;  bx[16] = 9;
    ax[14] = 4;  bx[14] = 10;
    ax[10] = 5;  bx[10] = 11;
    ax[12] = 6;  bx[12] = 12;

    ax[25] = 13; bx[7] = 17;
    ax[23] = 14; bx[5] = 18;
    ax[19] = 15; bx[1] = 19;
    ax[21] = 16; bx[3] = 20;

    ax[17] = 21; bx[17] = 25;
    ax[11] = 22; bx[11] = 26;
    ax[9] = 23;  bx[9] = 27;
    ax[15] = 24; bx[15] = 28;

    ax[26] = 33; bx[8] = 37;
    ax[20] = 34; bx[2] = 38;
    ax[18] = 35; bx[0] = 39;
    ax[24] = 36; bx[6] = 40;

    // Shared nodes.
    ax[4]    =   bx[22] = 7;
    
    ax[7]    =   bx[25] = 29;
    ax[5]    =   bx[23] = 30;
    ax[1]    =   bx[19] = 31;
    ax[3]    =   bx[21] = 32;

    ax[8]    =   bx[26] = 41;
    ax[2]    =   bx[20] = 42;
    ax[0]    =   bx[18] = 43;
    ax[6]    =   bx[24] = 44;

    // Initialize eltrans with these values.
/////    memcpy( eltrans.get_m_ctrl_idx().get_host_ptr(), ctrl_idx, 54*sizeof(int) );
/////    memcpy( eltrans.get_m_values().get_host_ptr(), grid_vals, 45*sizeof(float) );

    //-- Test ElTransQuery --//
/////    dray::ElTransQuery<dray::ElTrans_BernsteinShape<float,1,3>> eltransq;
/////    eltransq.m_eltrans = eltrans;
/////    eltransq.resize(num_queries);
/////    eltransq.m_el_ids = el_ids;
/////    eltransq.m_ref_pts = ref_pts;
/////    eltransq.query(active_idx);
/////    std::cout << "Test ElTransQuery" << std::endl;
/////    std::cout << "m_result_val";    eltransq.m_result_val.summary();
/////    std::cout << "m_result_deriv" << std::endl;
/////         eltransq.m_result_deriv.summary();
/////    std::cout << std::endl;

/////    //-- Test ElTransQuery2 --//
/////    typedef dray::ElTransQuery2<dray::ElTrans_BernsteinShape<float,1,3>,
/////                        dray::ElTrans_BernsteinShape<float,1,3>> ElTQ2;
/////    ElTQ2 eltransq2;
/////    eltransq2.m_q1.m_eltrans = eltrans;
/////    eltransq2.m_q2.m_eltrans = eltrans;
/////    eltransq2.resize(num_queries);
/////    eltransq2.m_q1.m_el_ids = el_ids;
/////    eltransq2.m_q2.m_el_ids = el_ids;
/////    eltransq2.m_q1.m_ref_pts = ref_pts;
/////    eltransq2.m_q2.m_ref_pts = ref_pts;
/////    eltransq2.query(active_idx);
/////    std::cout << "Test ElTransQuery2" << std::endl;
/////    std::cout << "m_result_val";    eltransq2.m_q1.m_result_val.summary(); eltransq2.m_q2.m_result_val.summary();
/////    std::cout << "m_result_deriv" << std::endl;
/////         eltransq2.m_q1.m_result_deriv.summary();
/////         eltransq2.m_q2.m_result_deriv.summary();
/////    std::cout << std::endl;
/////
/////    {
/////    ElTQ2::ptr_bundle_const_t ptrb = eltransq2.get_val_device_ptr_const();
/////    std::cout << "Some virtual values:  ";
/////    for (int ii = 0; ii < num_queries; ii++)
/////      std::cout << ElTQ2::get_val(ptrb, ii) << " ";
/////    std::cout << std::endl;
/////    }
/////
/////    {
/////    ElTQ2::ptr_bundle_const_t ptrb = eltransq2.get_deriv_device_ptr_const();
/////    std::cout << "Some virtual derivatives:  " << std::endl;
/////    for (int ii = 0; ii < num_queries; ii++)
/////      std::cout << ElTQ2::get_deriv(ptrb, ii);
/////    std::cout << std::endl;
/////    }


}
