#include "gtest/gtest.h"
#include "test_config.h"
#include <dray/camera.hpp>
#include <dray/triangle_mesh.hpp>
#include <dray/io/obj_reader.hpp>
#include <dray/utils/ray_utils.hpp>

TEST(dray_test, dray_test_unit)
{
  std::string file_name = std::string(DATA_DIR) + "unit_cube.obj";
  std::cout<<"File name "<<file_name<<"\n";
  
  dray::Array<dray::float32> vertices;
  dray::Array<dray::int32> indices;

  read_obj(file_name, vertices, indices);

  dray::TriangleMesh mesh(vertices, indices);
  dray::Camera camera;
  dray::Vec3f pos = dray::make_vec3f(10,10,10);
  dray::Vec3f look_at = dray::make_vec3f(5,5,5);
  camera.set_look_at(look_at);
  camera.set_pos(pos);
  camera.reset_to_bounds(mesh.get_bounds());
  dray::ray32 rays;
  camera.create_rays(rays);
  std::cout<<camera.print();
  mesh.intersect(rays);
  
  dray::save_depth(rays, camera.get_width(), camera.get_height());

}

//TEST(dray_test, dray_test_conference)
//{
//  std::string file_name = std::string(DATA_DIR) + "conference.obj";
//  std::cout<<"File name "<<file_name<<"\n";
//  
//  dray::Array<dray::float32> vertices;
//  dray::Array<dray::int32> indices;
//
//  read_obj(file_name, vertices, indices);
//
//  dray::TriangleMesh mesh(vertices, indices);
//  dray::Camera camera;
//
//  dray::Vec3f pos = dray::make_vec3f(30,19,5);
//  dray::Vec3f look_at = dray::make_vec3f(0,0,0);
//  dray::Vec3f up = dray::make_vec3f(0,0,1);
//
//  camera.set_look_at(look_at);
//  camera.set_pos(pos);
//  camera.set_up(up);
//  //camera.reset_to_bounds(mesh.get_bounds());
//  dray::ray32 rays;
//  camera.create_rays(rays);
//  std::cout<<camera.print();
//  mesh.intersect(rays);
// 
//  dray::save_depth(rays, camera.get_width(), camera.get_height());
//
//}
