#include "gtest/gtest.h"
#include "test_config.h"

#include "t_utils.hpp"
#include <dray/mfem2dray.hpp>
#include <dray/shaders.hpp>
#include <mfem.hpp>
#include <mfem/fem/conduitdatacollection.hpp>

#include <dray/camera.hpp>
#include <dray/utils/png_encoder.hpp>
#include <dray/utils/ray_utils.hpp>

#include <dray/math.hpp>

#include <fstream>
#include <stdlib.h>


// Helper function prototype.

// Returns pointer to new mesh and grid function.
// Caller is responsible to delete mesh_ptr and sol.
void construct_example_data(const int num_el,
                            mfem::Mesh *&mesh_ptr,
                            mfem::GridFunction * &sol,
                            int order = 2);

//
// TEST()
//
TEST(dray_volume_render, dray_volume_render_simple)
{
  std::string file_name = std::string(DATA_DIR) + "impeller/impeller";
  std::string output_path = prepare_output_dir();
  std::string output_file = conduit::utils::join_file_path(output_path, "impeller_vr");
  remove_test_image(output_file);

  mfem::Mesh *mfem_mesh_ptr;
  mfem::GridFunction *mfem_sol_ptr;

  // Initialize mfem data.
  //construct_example_data(50000, mfem_mesh_ptr, mfem_sol_ptr);
  //construct_example_data(1, mfem_mesh_ptr, mfem_sol_ptr);
  //mfem::ConduitDataCollection dcol("crazy_hex", mfem_mesh_ptr);
  //dcol.RegisterField("bananas", mfem_sol_ptr);
  //dcol.SetProtocol("conduit_bin");
  //dcol.SetCycle(0);
  //dcol.SetTime(0.0);
  //dcol.Save();

  mfem::ConduitDataCollection dcol(file_name);
  dcol.SetProtocol("conduit_bin");
  dcol.Load();
  mfem_mesh_ptr = dcol.GetMesh();
  mfem_sol_ptr = dcol.GetField("bananas");

  if (mfem_mesh_ptr->NURBSext)
  {
     mfem_mesh_ptr->SetCurvature(2);
  }
  mfem_mesh_ptr->GetNodes();

  // --- DRAY code --- //

  int space_P;
  dray::ElTransData<float,3> space_data = dray::import_mesh<float>(*mfem_mesh_ptr, space_P);

  int field_P;
  dray::ElTransData<float,1> field_data = dray::import_grid_function<float,1>(*mfem_sol_ptr, field_P);

  std::cout << "field_data.m_ctrl_idx ...   ";
  field_data.m_ctrl_idx.summary();
  std::cout << "field_data.m_values ...     ";
  field_data.m_values.summary();

  dray::Mesh<float> mesh(space_data, space_P);
  dray::Field<float> field(field_data, field_P);
  dray::MeshField<float> mesh_field(mesh, field);

  dray::ColorTable color_table("Spectral");
  color_table.add_alpha(0.f,  0.01f);
  color_table.add_alpha(0.1f, 0.09f);
  color_table.add_alpha(0.2f, 0.01f);
  color_table.add_alpha(0.3f, 0.09f);
  color_table.add_alpha(0.4f, 0.01f);
  color_table.add_alpha(0.5f, 0.01f);
  color_table.add_alpha(0.6f, 0.01f);
  color_table.add_alpha(0.7f, 0.09f);
  color_table.add_alpha(0.8f, 0.01f);
  color_table.add_alpha(0.9f, 0.01f);
  color_table.add_alpha(1.0f, 0.0f);
  dray::Shader::set_color_table(color_table);

  // Camera
  const int c_width = 1024;
  const int c_height = 1024;
  dray::Camera camera;
  camera.set_width(c_width);
  camera.set_height(c_height);
  camera.reset_to_bounds(mesh_field.get_bounds());
  dray::Array<dray::ray32> rays;
  camera.create_rays(rays);

  //
  // Volume rendering
  //

  float sample_dist;
  {
    constexpr int num_samples = 100;
    dray::AABB<> bounds = mesh_field.get_bounds();
    dray::float32 mag = (bounds.max() - bounds.min()).magnitude();
    sample_dist = mag / dray::float32(num_samples);
  }

  dray::Array<dray::Vec<dray::float32,4>> color_buffer = mesh_field.integrate(rays, sample_dist);

  {
    dray::PNGEncoder png_encoder;
    png_encoder.encode( (float *) color_buffer.get_host_ptr(), camera.get_width(), camera.get_height() );
    png_encoder.save(output_file + ".png");
    EXPECT_TRUE(check_test_image(output_file));
  }

#if 0
   //
   // Isosurface
   //

  camera.create_rays(rays);

  // Output isosurface, colorized by field spatial gradient magnitude.
  {
    float isovalues[5] = { 0.07, 0.005, 0, -8, -15 };
    const char* filenames[5] = {"isosurface_001.png",
                                "isosurface_+08.png",
                                "isosurface__00.png",
                                "isosurface_-08.png",
                                "isosurface_-15.png"};

    for (int iso_idx = 0; iso_idx < 1; iso_idx++)
    {
      std::cout<<"doing iso_surface "<<iso_idx<<" size "<<rays.size()<<"\n";
      dray::Array<dray::Vec4f> color_buffer = mesh_field.isosurface_gradient(rays, isovalues[iso_idx]);
      std::cout<<"done doing iso_surface "<<"\n";
      dray::PNGEncoder png_encoder;
      png_encoder.encode( (float *) color_buffer.get_host_ptr(), camera.get_width(), camera.get_height() );
      png_encoder.save(filenames[iso_idx]);

      printf("Finished rendering isosurface idx %d\n", iso_idx);
    }
  }
  // --- end DRAY  --- //
#endif
}


// --- MFEM code --- //

void construct_example_data(const int in_max_els,
                            mfem::Mesh *&out_mesh_ptr,
                            mfem::GridFunction * &out_sol_ptr,
                            int order)
{
  using namespace mfem;

  std::string file_name = std::string(DATA_DIR) + "beam-hex.mesh";
  //std::string file_name = std::string(DATA_DIR) + "beam-hex-nurbs.mesh";
  //std::string file_name = std::string(DATA_DIR) + "spiral_hex_p20.mesh";
  std::cout<<"File name "<<file_name<<"\n";

  Mesh *mesh = new Mesh(file_name.c_str(), 1, 1);
  int dim = mesh->Dimension();
  bool static_cond = false;
  int sdim = mesh->SpaceDimension();
  std::cout<<"Dim : "<<dim<<"\n"; //  Dims in referene space
  std::cout<<"Space Dim : "<<sdim<<"\n";

  const float max_els = in_max_els;
  // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   {
      int ref_levels =
         (int)floor(log(max_els/mesh->GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   mesh->ReorientTetMesh();

   // 4. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
   }
   else if (mesh->GetNodes())
   {
      fec = mesh->GetNodes()->OwnFEC();
      cout << "Using isoparametric FEs: " << fec->Name() << endl;
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
   }
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace->GetTrueVSize() << endl;

   // 5. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   Array<int> ess_tdof_list;
   if (mesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 6. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   LinearForm *b = new LinearForm(fespace);
   ConstantCoefficient one(1.0);
   b->AddDomainIntegrator(new DomainLFIntegrator(one));
   b->Assemble();

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction *_x = new GridFunction(fespace);
   GridFunction &x = *_x;
   x = 0.0;

   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));

   // 9. Assemble the bilinear form and the corresponding linear system,
   //    applying any necessary transformations such as: eliminating boundary
   //    conditions, applying conforming constraints for non-conforming AMR,
   //    static condensation, etc.
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   cout << "Size of linear system: " << A.Height() << endl;

   GSSmoother M(A);
   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);

   // 11. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);

   // 12. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
   sol_ofs.precision(8);
   x.Save(sol_ofs);


   // Output to arguments.
   out_mesh_ptr = mesh;
   out_sol_ptr = _x;

   printf("In construct_example(): fespace == %x\n", fespace);

   // TODO didn't there used to be some "delete" statements?
}
